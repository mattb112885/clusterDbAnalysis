#!/bin/bash

# Note - we only dump documentation for the scripts we expect people to use (i.e. not things in src/internal/)

mkdir "doc" 2> /dev/null
rm doc/help_texts

echo "o----------------------------------------------------------------------------------------------------------o" >> doc/help_texts
echo "|                                       COMMAND-LINE FUNCTIONS                                             |" >> doc/help_texts
echo "o----------------------------------------------------------------------------------------------------------o" >> doc/help_texts
echo "" >> doc/help_texts

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

cd utilities;

for file in *.py; do
    echo "*--------------------------------------------------------------------------------------------------------------------------*" >> ../../doc/help_texts
    echo "" >> ../../doc/help_texts
    echo "${file}" >> ../../doc/help_texts
    echo "" >> ../../doc/help_texts
    ${file} -h >> ../../doc/help_texts
done

cd ..;
cd ..;

echo "" >> doc/help_texts
echo "o----------------------------------------------------------------------------------------------------------o" >> doc/help_texts
echo "|                                       CONVENIENCE SCRIPTS                                                |" >> doc/help_texts
echo "o----------------------------------------------------------------------------------------------------------o" >> doc/help_texts
echo "" >> doc/help_texts

cd scripts;

for file in *.sh; do
    echo "*--------------------------------------------------------------------------------------------------------------------------*" >> ../doc/help_texts
    echo "" >> ../doc/help_texts
    echo "${file}" >> ../doc/help_texts
    echo "" >> ../doc/help_texts
    sh ${file} >> ../doc/help_texts
done