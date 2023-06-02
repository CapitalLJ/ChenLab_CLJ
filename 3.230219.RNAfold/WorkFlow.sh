mkdir fold ps pdf
for SUFFIX in {1..9}; do
  site0=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./mut.lst.tsv | cut -f 2)
  site1=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./mut.lst.tsv | cut -f 5)
  site2=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./mut.lst.tsv | cut -f 6)
  site3=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./mut.lst.tsv | cut -f 7)
  site4=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./mut.lst.tsv | cut -f 8)
  sed -n "${SUFFIX}p;${SUFFIX}q" ./mut.lst.tsv | perl ./snp4fasta.pl -f rawseq.fa --stdin >seq"${SUFFIX}".fa.tmp
  RNAfold -T 30 <seq"${SUFFIX}".fa.tmp >fold/seq"${SUFFIX}".fold
  RNAplot --pre "${site0} ${site0} 14 1 0 0 omark ${site1} ${site2} 8 0.33 1 0.33 omark ${site3} ${site4} 8 0 1 1 omark" <fold/seq"${SUFFIX}".fold
  rm seq"${SUFFIX}".fa.tmp
done

# shellcheck disable=SC2011
ls -- *_ss.ps | xargs -i{} mv {} ps/"${SUFFIX}"_{}
basename -s _ss.ps -a ps/*_ss.ps | xargs -i{} ps2pdf ps/{}_ss.ps pdf/{}.pdf
