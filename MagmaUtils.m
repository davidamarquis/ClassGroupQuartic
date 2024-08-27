RowspaceSupport := function(mat)
/* mat must be a SparseMatrix.BasisMatrix
I have assumed that BasisMatrix of Rowspace will be in HNF */
B:=BasisMatrix(Rowspace(mat));
return Support(B);
end function;

Diagonal := function(mat)
ncols:=NumberOfColumns(mat);
return [mat[i,i]: i in [1..ncols]];
end function;

ProbablePivotCols:=function(mat, pr)
  nrows := NumberOfRows(mat);
  ncols := NumberOfColumns(mat);
  mod_mats := KMatrixSpace(GF(pr), nrows, ncols);
  dense_mat:=Matrix(mat);
  mod_mat := mod_mats!dense_mat;
  ech_form:=EchelonForm(mod_mat);

  pivs:={};
  for rr in [1..nrows] do
    for cc in [rr..ncols] do
      if ech_form[rr,cc] ne 0 then
        pivs:=pivs join {cc};
        break;
      end if;
    end for;
  end for;
  return pivs;
end function;