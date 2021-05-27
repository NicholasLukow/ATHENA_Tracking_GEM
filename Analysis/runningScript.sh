#/bin/bash
while read file
do
	root -l -q './RunAnalyze.C("./'$file'")'
done < ./TreeList.txt
