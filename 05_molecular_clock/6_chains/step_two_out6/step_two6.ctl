          seed = -1
       seqfile = ../inp_data/metazoa_aln.phy
      treefile = ../inp_data/metazoa_tree_calib_MCMCtree.tree
      mcmcfile = mcmc6.txt
       outfile = out6.txt

         ndata = 1
       seqtype = 2
       usedata = 2 ../step_one_out/in.BV
         clock = 2
	 model = 3
	 alpha = 0.5 
         ncatG = 4

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1    * birth, death, sampling

   rgene_gamma = 2 5.1
  sigma2_gamma = 1 10

	 print = 1
        burnin = 100000
      sampfreq = 1000
       nsample = 20000
