mogrify -format pdf out/strain[0-9][0-9]_*.png 
pdftk out/strain[0-9][0-9]_*.pdf output out/strains.pdf
rm out/strain[0-9][0-9]_*.pdf
