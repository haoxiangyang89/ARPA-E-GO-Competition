include("def.jl");
include("dProc.jl");
include("mod_comp.jl");

fileAdd = "./Phase_0_Modified_RTS96/scenario_1/powersystem.raw";
genFile = "./Phase_0_Modified_RTS96/scenario_1/generator.csv";
contFile = "./Phase_0_Modified_RTS96/scenario_1/contingency.csv";

baseMVA,busSeg,loadSeg,shuntSeg,genSeg,branchSeg,transformerSeg = preProc(fileAdd);
busList,busDList,genList,genDList = busgenProc(baseMVA,busSeg,shuntSeg,loadSeg,genSeg,genFile);
brList,brListSingle,brDList = branchProc(baseMVA,branchSeg,transformerSeg);
contList,contDList = contProc(contFile);

fData = fixedData(baseMVA,busList,busDList,genList,genDList,brList,brListSingle,brDList);
uData = uncertainData(contList,contDList);

mp = buildMod(fData,uData);
solve(mp);
sphat = getvalue(mp[:sp0]);
sqhat = getvalue(mp[:sq0]);
spshat = getvalue(mp[:sp]);
sqshat = getvalue(mp[:sq]);
costhat = getobjectivevalue(mp);

open("solution1.txt","w") do f
  write(f, "--generation dispatch \n bus id,unit id,pg(MW),qg(MVar) \n");
  for i in fData.genList
    loc = fData.genDList[i].Loc;
    name = fData.genDList[i].Name;
    spTemp = sphat[i]*fData.baseMVA;
    sqTemp = sqhat[i]*fData.baseMVA;
    write(f, "$loc,$name,$spTemp,$sqTemp \n");
  end
  write(f,"--end of generation dispatch \n");
end

open("solution2.txt","w") do f
  write(f, "--contingency generator \n conID,genID,busID,unitID,q(MW) \n");
  for s in uData.contList
    for i in fData.genList
      loc = fData.genDList[i].Loc;
      idTemp = fData.genDList[i].ID;
      name = fData.genDList[i].Name;
      sqTemp = sqshat[i,s]*fData.baseMVA;
      write(f,"$s,$idTemp,$loc,$name,$sqTemp \n");
    end
  end
  write(f,"--end of contingency generator \n --bus \n contingency id,bus id,v(pu),theta(deg) \n");
  for i in fData.busList
    id = fData.busDList[i].ID;
    name = fData.busDList[i].Name;
    vTemp = getvalue(mp[:v0][i]);
    θTemp = getvalue(mp[:v0][i]/pi*180);
    write(f,"0,$id,$name,$vTemp,$θTemp \n");
  end
  for s in uData.contList
    for i in fData.busList
      id = fData.busDList[i].ID;
      name = fData.busDList[i].Name;
      vTemp = getvalue(mp[:v][i,s]);
      θTemp = getvalue(mp[:v][i,s]/pi*180);
      write(f,"$s,$id,$name,$vTemp,$θTemp \n");
    end
  end
  write(f,"--end of bus \n --Delta \n contingency id,Delta(MW) \n");
  for s in uData.contList
    pΔTemp = getvalue(mp[:pΔ][s])*fData.baseMVA;
    write(f,"$s,$pΔTemp \n");
  end
  write(f,"--end of Delta \n --line flow \n contingency id,line id,origin bus id,destination bus id,circuit id,p_origin(MW),q_origin(MVar),p_destination(MW),q_destination(MVar) \n");
  for i in fData.brListSingle
    idTemp = fData.brDList[i].ID;
    revidTemp = fData.brDList[i].revID;
    fromTemp = fData.brDList[i].From;
    toTemp = fData.brDList[i].To;
    name = fData.brDList[i].CKT;
    pTemp = getvalue(mp[:p0][i])*fData.baseMVA;
    qTemp = getvalue(mp[:q0][i])*fData.baseMVA;
    pRevTemp = getvalue(mp[:p0][revidTemp])*fData.baseMVA;
    qRevTemp = getvalue(mp[:q0][revidTemp])*fData.baseMVA;
    write(f,"0,$idTemp,$fromTemp,$toTemp,$name,$pTemp,$qTemp,$pRevTemp,$qRevTemp \n");
  end
  for s in uData.contList
    for i in fData.brListSingle
      idTemp = fData.brDList[i].ID;
      revidTemp = fData.brDList[i].revID;
      fromTemp = fData.brDList[i].From;
      toTemp = fData.brDList[i].To;
      name = fData.brDList[i].CKT;
      pTemp = getvalue(mp[:p][i,s])*fData.baseMVA;
      qTemp = getvalue(mp[:q][i,s])*fData.baseMVA;
      pRevTemp = getvalue(mp[:p][revidTemp,s])*fData.baseMVA;
      qRevTemp = getvalue(mp[:q][revidTemp,s])*fData.baseMVA;
      write(f,"$s,$idTemp,$fromTemp,$toTemp,$name,$pTemp,$qTemp,$pRevTemp,$qRevTemp \n");
    end
  end
  write(f,"--end of line flow \n")
end
