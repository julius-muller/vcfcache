# Sort records only within a CHROM+POS block. A positional annotator can reorder
# records that share a position (multi-allelic split order, indel vs SNV tie-break)
# but never move one across positions, so this yields an order-insensitive digest
# in one streaming pass instead of an external sort over millions of records.
function flush(  i,j,t){ for(i=2;i<=n;i++){t=b[i];for(j=i-1;j>=1&&b[j]>t;j--)b[j+1]=b[j];b[j+1]=t} for(i=1;i<=n;i++) print b[i]; n=0 }
{ p=index($0,"\t"); q=index(substr($0,p+1),"\t"); k=substr($0,1,p+q)
  if(k!=pk){ flush(); pk=k } b[++n]=$0 }
END{ flush() }
