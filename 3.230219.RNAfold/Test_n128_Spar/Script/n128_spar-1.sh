#!/usr/bin bash
for SUFFIX in {1..10000}; do
  site0=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./n128_Spar_SNP.tsv | cut -f 2)
  site1=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./n128_Spar_SNP.tsv | cut -f 1)
  site2=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./n128_Spar_SNP.tsv | cut -f 3)
  site3=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./n128_Spar_SNP.tsv | cut -f 4)
  # site1=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./n128_Spar_SNP.tsv | cut -f 5)
  # site2=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./n128_Spar_SNP.tsv | cut -f 6)
  # site3=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./n128_Spar_SNP.tsv | cut -f 7)
  # site4=$(sed -n "${SUFFIX}p;${SUFFIX}q" ./n128_Spar_SNP.tsv | cut -f 8)
  sed -n "${SUFFIX}p;${SUFFIX}q" ./n128_Spar_SNP.tsv | perl ./snp4fasta.pl -f n128_Spar.gene.fa --stdin >seq"${SUFFIX}".fa.tmp
  RNAfold -T 30 <seq"${SUFFIX}".fa.tmp >fold/"${site1}"_"${site0}"_"${site2}"_"${site3}".fold
  # RNAplot --pre "${site0} ${site0} 14 1 0 0 omark ${site1} ${site2} 8 0.33 1 0.33 omark ${site3} ${site4} 8 0 1 1 omark" <fold/seq"${SUFFIX}".fold
  RNAplot --pre "${site0} ${site0} 14 1 0 0 omark" <fold/"${site1}"_"${site0}"_"${site2}"_"${site3}".fold
  rm seq"${SUFFIX}".fa.tmp
  # shellcheck disable=SC2011

done
# ls -- *_ss.ps | xargs -i{} mv {} ps/{}
# # basename -s _ss.ps -a ps/*_ss.ps | xargs -i{} ps2pdf ps/{}_ss.ps pdf/{}.pdf
# ls ps/*_ss.ps | xargs -I {} sh -c 'ps2pdf {} pdf/$(basename {} _ss.ps).pdf'


