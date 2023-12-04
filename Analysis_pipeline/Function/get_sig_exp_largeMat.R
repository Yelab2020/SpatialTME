as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
get_sig_exp_largeMat <- function (se.obj = WholeTissueSC, DefineTypes = "Majortypes", 
                                  sig_scran = sig_scran) 
{
  mat <- as_matrix(se.obj@assays$RNA@data)
  norm_exp <- 2^mat - 1
  id <- se.obj@meta.data[, DefineTypes]
  ExprSubset <- norm_exp[sig_scran, ]
  Sig_exp <- NULL
  for (i in unique(id)) {
    Sig_exp <- cbind(Sig_exp, (apply(ExprSubset, 1, function(y) mean(y[which(id == 
                                                                               i)]))))
  }
  colnames(Sig_exp) <- unique(id)
  return(Sig_exp)
}