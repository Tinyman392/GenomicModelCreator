#!/bin/bash
# convert_beta_new.sh [beta model dir] <paired end>

if [ -d .stage ]; then
	rm -rf .stage
fi

mDir=$(echo $1 | sed 's/\/$//g')
mkdir .stage

pEnd=False
if [ $# -gt 2 ]; then
	pEnd=$3
fi

cp $mDir/attrOrder .stage/model.attrOrder
> .stage/model.labels.map
if [ -f $mDir/classes_map.tab ]; then
	cp $mDir/classes_map.tab .stage/model.labels.map
# else
# 	for i in $(cut -f3 -d$'\t' $mDir/model*0.true | sort -u); do
# 		echo -e $i'\t'$i
# 	done > .stage/model.labels.map
fi
mkdir .stage/all
cp $mDir/* .stage/all/

kLen=$(expr length "$(head -1 $mDir/attrOrder | cut -f1 -d$'\t')")
prAb=False
if [ $kLen -gt 12 ]; then
	prAb=True
fi

s='{"pairedEnd": '$pEnd', "kmerSize": '$kLen', "presence_absence": '$prAb'}'
echo $s > .stage/model.params

cp -r .stage/* $mDir/ 