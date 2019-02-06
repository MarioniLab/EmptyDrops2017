# Get the data
wget https://jmlab-gitlab.cruk.cam.ac.uk/publications/EmptyDrops2017-DataFiles/raw/08a775d772aadd25404f141bf95fce79412a0f3b/placenta_vento-tormo.zip
unzip placenta_vento-tormo.zip

for x in placenta1 placenta2 placenta3 placenta4 placenta5 placenta6
do
    mkdir -p $x
    tar -xvf $x"_raw_gene_bc_matrices.tar.gz" -C $x
    rm $x"_raw_gene_bc_matrices.tar.gz"
done
rm placenta_vento-tormo.zip

# Download the fetal/maternal annotation.
wget https://jmlab-gitlab.cruk.cam.ac.uk/publications/EmptyDrops2017-DataFiles/raw/master/cell_origin.csv
