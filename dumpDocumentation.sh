#!/bin/bash

mkdir "doc" 2> /dev/null
rm doc/help_texts

cd src;
for file in *.py; do
    echo "*--------------------------------------------------------------------------------------------------------------------------*" >> ../doc/help_texts
    echo "" >> ../doc/help_texts
    echo "${file}" >> ../doc/help_texts
    echo "" >> ../doc/help_texts
    ${file} -h >> ../doc/help_texts
done

for file in *.sh; do
    echo "*--------------------------------------------------------------------------------------------------------------------------*" >> ../doc/help_texts
    echo "" >> ../doc/help_texts
    echo "${file}" >> ../doc/help_texts
    echo "" >> ../doc/help_texts
    sh ${file} >> ../doc/help_texts
done

cd ..;

cd scripts;

for file in *.sh; do
    echo "*--------------------------------------------------------------------------------------------------------------------------*" >> ../doc/help_texts
    echo "" >> ../doc/help_texts
    echo "${file}" >> ../doc/help_texts
    echo "" >> ../doc/help_texts
    sh ${file} >> ../doc/help_texts
done