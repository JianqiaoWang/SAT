Find.plink.block =  function(snplist, ref.bed, output.dir = "./temp/"){
# find the block based on the predefined independent LD block
  dir.create(output.dir)
  write.csv(snplist ,row.names = F, quote = FALSE, file = paste(output.dir, "target.snplist", sep = ""))

  Block.map = data.table::fread("../supp/ld_block_hg19.txt")
  snp.pos =  data.table::fread("../supp/sumstat.pos.b37.txt") %>%
    dplyr::slice(match(snplist, variant_id)) %>%
    dplyr::arrange(CHR, POS) %>%
    dplyr::mutate(CHR = paste0("chr",CHR))
  snp.pos$block = NA

  for(i in 1:nrow(snp.pos)){
    snp.pos$block[i] = which(snp.pos$CHR[i] == Block.map$chr & (snp.pos$POS[i] >= Block.map$start) & (snp.pos$POS[i] < Block.map$stop))
  }

  Block.list = split(snp.pos$variant_id, snp.pos$block)

  return(Block.list)
}

Find.plink.block.elaborate =  function(snplist, ref.bed, output.dir = "./temp/"){
  # find the block based on the predefined independent LD block
  dir.create(output.dir)
  write.csv(snplist ,row.names = F, quote = FALSE, file = paste(output.dir, "target.snplist", sep = ""))

  Block.map = data.table::fread("../supp/ld_block_hg19.txt")
  snp.pos =  data.table::fread("../supp/sumstat.pos.b37.txt") %>%
    dplyr::slice(match(snplist, variant_id)) %>%
    dplyr::arrange(CHR, POS) %>%
    dplyr::mutate(CHR = paste0("chr",CHR))
  snp.pos$block = NA

  for(i in 1:nrow(snp.pos)){
    snp.pos$block[i] = which(snp.pos$CHR[i] == Block.map$chr & (snp.pos$POS[i] >= Block.map$start) & (snp.pos$POS[i] < Block.map$stop))
  }

  block.map1 =  snp.pos %>% group_by(block, CHR) %>%
    dplyr::summarise(num = n(), start = min(POS), end = max(POS) )

  block.map1$flag = 1

  for(i in 2:nrow(block.map1)){

    if( (block.map1$start[i] - block.map1$end[i-1]) < 100000 & (block.map1$CHR[i] == block.map1$CHR[i-1])  ){
      block.map1$flag[i] =  block.map1$flag[i-1]
    }else{
      block.map1$flag[i] =  block.map1$flag[i-1] + 1
    }
  }
  Map1 = setNames(block.map1$flag, as.character(block.map1$block))
  snp.pos$block2 = Map1[as.character(snp.pos$block)]
  Block.list = split(snp.pos$variant_id, snp.pos$block2)
  return(Block.list)
}

Find.plink.block.loose =  function(snplist, ref.bed, output.dir = "./temp/"){

  dir.create(output.dir)
  write.csv(snplist ,row.names = F, quote = FALSE, file = paste(output.dir, "target.snplist", sep = ""))

  txt = paste("plink --bfile ", ref.bed,  " --extract ",
              paste(output.dir, "target.snplist", sep = ""),
              " --blocks no-pheno-req no-small-max-span --blocks-max-kb 2000 --blocks-min-maf 0.005 --blocks-strong-lowci 0.5005 --blocks-strong-highci 0.8305 --blocks-recomb-highci 0.80 --out ", output.dir,"target.block", sep = "")

  system(paste("module load plink \n", txt, sep = ""))

  Block = fread(paste(output.dir,"target.block.blocks.det", sep = ""))
  Block.list = strsplit(Block$SNPS, split = "|", fixed = T)

  if(sum(Block$NSNPS) == length(snplist)){ print("All variants are classified into plink estimated haplotype blocks") }

  if(sum(Block$NSNPS) != length(snplist)){

    print("Some variants are not included in the blocks")
    loc.start.end = sapply( Block.list, function(x){
      loc.start.end.info = match( c(head(x, 1), tail(x, 1)  ), snplist)
      return( c(  min(loc.start.end.info), max(loc.start.end.info)  ))
    })

    # from close interval [a1, a2], [b1, b2] , convert it to [1,a1) [a1, b1), [b1, tail]
    loc.vector = as.vector(t(loc.start.end)[,1]) # vectorize by row
    complete.block = data.frame(start = c(1, loc.vector), end =  c(loc.vector, length(snplist)+1))
    complete.block = complete.block[! (complete.block$start >= complete.block$end), ]
    Block.list = Map(function(i,j) snplist[i:j], complete.block$start  , complete.block$end-1 )

    if((length(unlist(Block.list)) == length(snplist))  ){
      print("All variants are divided into blocks")
    }
  }

  return(Block.list)

}


Find.plink.block.original =  function(snplist, ref.bed, output.dir = "./temp/"){
  # too conservative
  dir.create(output.dir)
  write.csv(snplist ,row.names = F, quote = FALSE, file = paste(output.dir, "target.snplist", sep = ""))

  txt = paste("plink --bfile ", ref.bed,  " --extract ",
              paste(output.dir, "target.snplist", sep = ""),
              " --blocks no-pheno-req no-small-max-span --blocks-max-kb 2000 --blocks-min-maf 0.005 --blocks-strong-lowci 0.5005 --blocks-strong-highci 0.8305 --blocks-recomb-highci 0.80 --out ", output.dir,"target.block", sep = "")

  system(paste("module load plink \n", txt, sep = ""))

  Block = fread(paste(output.dir,"target.block.blocks.det", sep = ""))
  Block.list = strsplit(Block$SNPS, split = "|", fixed = T)

  if(sum(Block$NSNPS) == length(snplist)){ print("All variants are classified into plink estimated haplotype blocks") }

  if(sum(Block$NSNPS) != length(snplist)){

    print("Some variants are not included in the blocks")
    loc.start.end = sapply( Block.list, function(x){
      loc.start.end.info = match( c(head(x, 1), tail(x, 1)  ), snplist)
      return( c(  min(loc.start.end.info), max(loc.start.end.info)  ))
    })

    # from close interval [a, b] to left close right open [a, b+1)
    loc.start.end.right.open = t(loc.start.end) ;  loc.start.end.right.open[,2] = loc.start.end.right.open[,2]+1;
    loc.vector = as.vector(t(loc.start.end.right.open)) # vectorize by row
    complete.block = data.frame(start = c(1, loc.vector), end =  c(loc.vector, length(snplist)+1))
    complete.block = complete.block[! (complete.block$start >= complete.block$end), ]
    Block.list.sub = Block.list
    Block.list = Map(function(i,j) snplist[i:j], complete.block$start  , complete.block$end-1 )

    if(all(Block.list.sub %in% Block.list) & (length(unlist(Block.list)) == length(snplist))  ){
      print("All variants are divided into blocks")
    }
  }

  return(Block.list)

}

#
# GenerateBlock.plink = function(snplist, ref.geno, output.dir){
#
#     write.csv(snplist ,row.names = F, quote = FALSE, file = paste(output.dir, "target.snplist", sep = ""))
#
#     txt = paste("plink --bfile ", ref.geno,  " --extract ",
#     paste(output.dir, "target.snplist", sep = ""),
#     " --blocks no-pheno-req no-small-max-span --blocks-max-kb 2000 --blocks-min-maf 0.001 --blocks-strong-lowci 0.5005 --blocks-strong-highci 0.8305 --blocks-recomb-highci 0.80 --out ", output.dir,"target.block", sep = "")
#
#     a <- system(paste("module load plink \n", txt, sep = ""))
#
#     Block = fread(paste(output.dir,"target.block.blocks.det", sep = ""))
#
#     Block.list = strsplit(Block$SNPS, split = "|", fixed = T)
#
#     return(Block.list)
# }

