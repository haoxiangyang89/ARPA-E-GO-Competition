# cost calculation functions
function costCal(fData,sphat,sqhat)
  cost = sum(sum(fData.genDList[l].cParams[n]*(sphat[l]*fData.baseMVA)^n for n in fData.genDList[l].cn) for l in fData.genList);
  return cost
end
