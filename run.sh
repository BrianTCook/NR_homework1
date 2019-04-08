#!/bin/bash

echo ""
echo "-------------------------------------"
echo "Brian Cook's homework 1 solutions"
echo "-------------------------------------"

echo ""
echo "-------------------------------------"
echo "Running the Python scripts for each part"
echo "-------------------------------------"

echo ""
echo "-------------------------------------"
echo "Problem 1, Part A"
echo "-------------------------------------"

python3 homework1problem1parta.py

echo ""
echo "-------------------------------------"
echo "Problem 1, Part B"
echo "-------------------------------------"

python3 homework1problem1partb.py

echo ""
echo "-------------------------------------"
echo "Problem 2, Part A"
echo "-------------------------------------"

python3 homework1problem2parta.py

echo ""
echo "-------------------------------------"
echo "Problem 2, Part B"
echo "-------------------------------------"

python3 homework1problem2partb.py

echo ""
echo "-------------------------------------"
echo "Problem 2, Part C"
echo "-------------------------------------"

python3 homework1problem2partc.py

echo ""
echo "-------------------------------------"
echo "Problem 2, Part D"
echo "-------------------------------------"

python3 homework1problem2partd.py

echo ""
echo "-------------------------------------"
echo "Problem 2, Part E"
echo "-------------------------------------"

python3 homework1problem2parte.py

echo ""
echo "-------------------------------------"
echo "Problem 2, Part F"
echo "-------------------------------------"

python3 homework1problem2partf.py

echo ""
echo "-------------------------------------"
echo "Problem 2, Part G"
echo "-------------------------------------"

python3 homework1problem2partg.py

echo ""
echo "-------------------------------------"
echo "Problem 2, Part H"
echo "-------------------------------------"

python3 homework1problem2parth.py

echo ""
echo "-------------------------------------"
echo "Downloading data for Problem 3"
echo "-------------------------------------"

wget https://home.strw.leidenuniv.nl/~daalen/files/satgals_m11.txt
wget https://home.strw.leidenuniv.nl/~daalen/files/satgals_m12.txt
wget https://home.strw.leidenuniv.nl/~daalen/files/satgals_m13.txt
wget https://home.strw.leidenuniv.nl/~daalen/files/satgals_m14.txt
wget https://home.strw.leidenuniv.nl/~daalen/files/satgals_m15.txt

echo ""
echo "-------------------------------------"
echo "Problem 3, Part A"
echo "-------------------------------------"

echo "Current Python script cannot apply the simplex method to find the optimal values of (a,b,c)"
echo "Would need to work some more with the dictionary {(a,b,c): neglogL(a,b,c)} to get this to work."
#python3 homework1problem3parta.py

echo ""
echo "-------------------------------------"
echo "Problem 3, Part B"
echo "-------------------------------------"

echo "I don't have a Python script for this part"

echo ""
echo "-------------------------------------"
echo "Generating the pdf"
echo "-------------------------------------"

pdflatex homework1_cook.tex
bibtex homework1_cook.aux
pdflatex homework1_cook.tex
pdflatex homework1_cook.tex

echo ""
echo "-------------------------------------"
echo "My code/outputs are saved as homework1_cook.pdf"
echo "-------------------------------------"
