# script to extract krakekn results across different taxonomix levels
# repuires lauch in working directory
# results in: 2-kraken2/sumKrakTest_*.txt

# Read input from CL
man_switch="off"
while [ $# != "0" ] ; do
 if [[ ${switch} == 'off' ]]; then
  input=$1
  man_switch="on"
 else
  input="$input $1"
 fi
shift
done


for f in ${input}; do
 sid=$(echo $f |sed -e 'sX2-kraken2/XX' -e 'sX_krakenRep.txtXX')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="S"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 2-kraken2/sumKrakTest_species.txt

for f in ${input}; do
 sid=$(echo $f |sed -e 'sX2-kraken2/XX' -e 'sX_krakenRep.txtXX')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="G"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 2-kraken2/sumKrakTest_genus.txt

for f in ${input}; do
 sid=$(echo $f |sed -e 'sX2-kraken2/XX' -e 'sX_krakenRep.txtXX')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="F"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 2-kraken2/sumKrakTest_family.txt

for f in ${input}; do
 sid=$(echo $f |sed -e 'sX2-kraken2/XX' -e 'sX_krakenRep.txtXX')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="O"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 2-kraken2/sumKrakTest_order.txt

for f in ${input}; do
 sid=$(echo $f |sed -e 'sX2-kraken2/XX' -e 'sX_krakenRep.txtXX')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="C"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 2-kraken2/sumKrakTest_class.txt

for f in ${input}; do
 sid=$(echo $f |sed -e 'sX2-kraken2/XX' -e 'sX_krakenRep.txtXX')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="P"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 2-kraken2/sumKrakTest_phylum.txt

for f in ${input}; do
 sid=$(echo $f |sed -e 'sX2-kraken2/XX' -e 'sX_krakenRep.txtXX')
 cat $f |awk '{$6=$6}1' |awk '{if($4=="D"){print $0}}' |cut -d' ' -f1,2,6- |sed 's/ /@/' |sed 's/ /_/g' |sed 's/@/ /' |sed 's/_/ /' |while read line; do echo $line $sid; done
done > 2-kraken2/sumKrakTest_domain.txt
