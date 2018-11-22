# get the data
wget https://ndownloader.figshare.com/articles/7373870/versions/1
unzip 1 # it comes zipped in this thing

for x in placenta1 placenta2 placenta3 placenta4 placenta5 placenta6
do
    mkdir -p $x
    tar -xvf $x"_raw_gene_bc_matrices.tar.gz" -C $x
    rm $x"_raw_gene_bc_matrices.tar.gz"
done
rm 1
