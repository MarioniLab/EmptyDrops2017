for x in \
    http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz \
    http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neurons_900/neurons_900_raw_gene_bc_matrices.tar.gz \
    http://cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/293t_raw_gene_bc_matrices.tar.gz \
    http://cf.10xgenomics.com/samples/cell-exp/1.1.0/jurkat/jurkat_raw_gene_bc_matrices.tar.gz \
    http://cf.10xgenomics.com/samples/cell-exp/2.1.0/t_4k/t_4k_raw_gene_bc_matrices.tar.gz \
    http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neuron_9k/neuron_9k_raw_gene_bc_matrices.tar.gz
do
    wget $x
    destname=$(basename $x) 
    stub=$(echo $destname | sed "s/_raw_.*//")
    mkdir -p $stub
    tar -xvf $destname -C $stub
    rm $destname
done
