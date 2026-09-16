"""chainmap.py -- minimal, from-scratch UCSC chain file parser + coordinate/interval mapper.

Written instead of relying solely on `pyliftover` because this pipeline needs BOTH directions:
assembly -> hg38 (the main CpG-projection direction) and hg38 -> assembly (to restrict a modbed
tabix query to one small region without parsing a whole haplotype). `pyliftover` only exposes
target -> query. Also lets point-mapping be vectorized later (np.searchsorted) if a per-CpG
Python loop turns out to be the Stage 2 bottleneck at whole-chromosome scale -- not done yet,
since it hasn't been measured as one (see PRIORITIES-style lesson from ont_asm_caller: measure
before optimizing).

Correctness is NOT assumed. `validate_against_pyliftover()` below cross-checks a sample of
mapped points against `pyliftover` (a mature, independent implementation) -- run this on every
new chain file before trusting its bulk output for anything.

Confirmed from a real chain header (HG00097_hap1_vs_GRCh38.chain.gz):
    chain 22412982147 CM094061.1 242754100 + 4781 242753557 chr2 242193529 + 10198 242183504 1
`target` (t*) = the assembly contig (bare accession, e.g. "CM094061.1" -- NOT the
"SAMPLE#hap#accession" naming modbed reads use, strip that prefix before lookup).
`query` (q*) = GRCh38/CHM13. So target -> query is assembly -> reference, the direction this
pipeline mostly needs.

UCSC chain format spec: https://genome.ucsc.edu/goldenPath/help/chain.html
    chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    size dt dq       (repeated)
    size             (last line of a block, no dt/dq)
    <blank line separates chains>

Coordinates are 0-based half-open, given relative to the sequence's OWN strand as named by
tStrand/qStrand -- if that strand is '-', the given start/end are relative to the reverse
complement, so the true plus-strand coordinate is (size - given). Normalized to always be
plus-strand-relative in `parse_chain()`'s output so callers never have to think about it again.
"""
import bisect
import gzip


def _open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def parse_chain(path):
    """Return a list of ungapped-block tuples:
    (t_name, t0, t1, q_name, q0, q1, flip) -- t0/t1/q0/q1 always plus-strand-relative;
    flip=True means the block is on opposite strands (t_strand != q_strand), so within the
    block, t and q positions run in opposite directions."""
    blocks = []
    header = None
    t_pos = q_pos = None
    with _open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                header = None
                continue
            if line.startswith('chain'):
                f = line.split()
                header = dict(t_name=f[2], t_size=int(f[3]), t_strand=f[4], t_start=int(f[5]),
                              q_name=f[7], q_size=int(f[8]), q_strand=f[9], q_start=int(f[10]))
                t_pos, q_pos = header['t_start'], header['q_start']
                continue
            if header is None:
                continue
            f = line.split()
            size = int(f[0])
            t0, q0 = t_pos, q_pos
            t1, q1 = t_pos + size, q_pos + size
            if header['t_strand'] == '-':
                nt0, nt1 = header['t_size'] - t1, header['t_size'] - t0
            else:
                nt0, nt1 = t0, t1
            if header['q_strand'] == '-':
                nq0, nq1 = header['q_size'] - q1, header['q_size'] - q0
            else:
                nq0, nq1 = q0, q1
            blocks.append((header['t_name'], nt0, nt1, header['q_name'], nq0, nq1,
                            header['t_strand'] != header['q_strand']))
            if len(f) == 3:
                t_pos += size + int(f[1])
                q_pos += size + int(f[2])
            else:
                t_pos += size
                q_pos += size
    return blocks


class ChainMap:
    def __init__(self, chain_path):
        self.by_t = {}
        self.by_q = {}
        for t_name, t0, t1, q_name, q0, q1, flip in parse_chain(chain_path):
            self.by_t.setdefault(t_name, []).append((t0, t1, q_name, q0, q1, flip))
            self.by_q.setdefault(q_name, []).append((q0, q1, t_name, t0, t1, flip))
        for d in (self.by_t, self.by_q):
            for k in d:
                d[k].sort()

    @staticmethod
    def _point(index, name, pos):
        recs = index.get(name)
        if not recs:
            return None
        i = bisect.bisect_right(recs, (pos,)) - 1
        # bisect on tuples starting with the same field we sorted by (s0) works since s0 is
        # the first tuple element in both by_t/by_q's stored records
        if i < 0:
            return None
        s0, s1, oname, o0, o1, flip = recs[i]
        if not (s0 <= pos < s1):
            return None
        offset = pos - s0
        return (oname, (o1 - 1 - offset) if flip else (o0 + offset))

    def t_to_q(self, t_name, pos):
        """assembly contig coordinate -> reference (GRCh38/CHM13) coordinate."""
        return self._point(self.by_t, t_name, pos)

    def q_to_t(self, q_name, pos):
        """reference coordinate -> assembly contig coordinate."""
        return self._point(self.by_q, q_name, pos)

    def q_interval_to_t(self, q_name, start, end):
        """All target-side sub-intervals whose query-side block overlaps [start, end).
        Used to find the assembly-coordinate region to tabix-query for a given hg38 region --
        NOT assumed to be a single contiguous block (indels/rearrangements can split it)."""
        out = []
        for q0, q1, t_name, t0, t1, flip in self.by_q.get(q_name, []):
            lo, hi = max(q0, start), min(q1, end)
            if lo < hi:
                if flip:
                    tlo, thi = t0 + (q1 - hi), t0 + (q1 - lo)
                else:
                    tlo, thi = t0 + (lo - q0), t0 + (hi - q0)
                out.append((t_name, tlo, thi))
        return out


def validate_against_pyliftover(chain_path, cm, t_name, points, n=20):
    """Cross-check `cm.t_to_q` against pyliftover's independent implementation on up to `n` of
    `points` (assembly-contig positions). Returns (n_checked, n_agree, disagreements)."""
    from pyliftover import LiftOver
    lo = LiftOver(chain_path)
    checked, agree, bad = 0, 0, []
    for pos in points[:n]:
        mine = cm.t_to_q(t_name, pos)
        theirs = lo.convert_coordinate(t_name, pos)
        checked += 1
        their_result = (theirs[0][0], theirs[0][1]) if theirs else None
        if mine == their_result:
            agree += 1
        else:
            bad.append((pos, mine, their_result))
    return checked, agree, bad
