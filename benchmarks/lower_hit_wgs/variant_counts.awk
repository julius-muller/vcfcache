BEGIN { FS = OFS = "\t" }
{
    kind = (length($1) == 1 && length($2) == 1 ? "SNP" : "indel")
    all_count++; total[kind]++
    if ($3 != "." && $3 + 0 >= 30) { high_count++; high[kind]++ }
}
END {
    print "all", "all", all_count + 0
    print "SNP", "all", total["SNP"] + 0
    print "indel", "all", total["indel"] + 0
    print "all", "QUAL>=30", high_count + 0
    print "SNP", "QUAL>=30", high["SNP"] + 0
    print "indel", "QUAL>=30", high["indel"] + 0
}
