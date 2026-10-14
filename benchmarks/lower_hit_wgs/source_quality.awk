# Coalesce high-confidence gVCF intervals while counting source call quality.
# Input: CHROM POS END REF ALT QUAL FILTER GT DP GQ MIN_DP.
BEGIN { FS = OFS = "\t"; records = 0; gqm = gql = dpm = dpl = qm = ql = 0 }
function flush() {
    if (active) print chromosome, start, finish > bed
}
{
    c = $1; sub(/^chr/, "", c)
    if (c !~ /^[0-9]+$/ || c < 1 || c > 22) next
    if ($7 != "." && $7 != "PASS") next
    gq = ($10 == "." ? -1 : $10 + 0)
    dp = ($11 != "." ? $11 + 0 : ($9 != "." ? $9 + 0 : -1))
    if ($8 ~ /[1-9]/) {
        records++
        gqm += ($10 == ".")
        gql += ($10 != "." && gq < 20)
        dpm += ($9 == ".")
        dpl += ($9 != "." && $9 + 0 < 10)
        qm += ($6 == ".")
        ql += ($6 != "." && $6 + 0 < 30)
    }
    if (gq >= 20 && dp >= 10 && $8 !~ /\./) {
        current_chromosome = "chr" c
        lo = $2 - 1
        hi = ($3 == "." ? lo + length($4) : $3 + 0)
        if (active && current_chromosome == chromosome && lo <= finish) {
            if (hi > finish) finish = hi
        } else {
            flush()
            chromosome = current_chromosome; start = lo; finish = hi; active = 1
        }
    }
}
END { flush(); print records, gqm, gql, dpm, dpl, qm, ql }
