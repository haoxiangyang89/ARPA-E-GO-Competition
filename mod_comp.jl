# model of the ARPA-E Competition
using Ipopt, Complementarity;

function buildMod(fData,uData)
  # Input:
  # fData - grid data with all parameters
  # uData - contingency data

  # readin the data from the structure
  baseMVA = fData.baseMVA;
  bList = fData.busList;
  bData = fData.busDList;

  gList = fData.genList;
  gData = fData.genDList;

  brList = fData.brList;
  brData = fData.brDList;

  S = uData.contList;
  contData = uData.contDList;

  # set up the model
  mp = Model(solver = IpoptSolver());

  # create the variables for the base case
  @variable(mp,bData[i].Vmin <= v0[i in bList] <= bData[i].Vmax);
  @variable(mp,gData[l].Pmin <= sp0[l in gList] <= gData[l].Pmax);
  @variable(mp,gData[l].Qmin <= sq0[l in gList] <= gData[l].Qmax);
  @variable(mp,p0[k in brList]);
  @variable(mp,q0[k in brList]);
  @variable(mp,psh0[i in bList]);
  @variable(mp,qsh0[i in bList]);
  @variable(mp,θ0[i in bList]);

  # create the constraints for the base case
  @NLconstraint(mp,flowBound0[k in brList],p0[k]^2 + q0[k]^2 <= brData[k].t);
  @NLconstraint(mp,pShunt0[i in bList],psh0[i] == bData[i].gsh*v0[i]^2);
  @NLconstraint(mp,qShunt0[i in bList],qsh0[i] == -bData[i].bsh*v0[i]^2);
  @NLconstraint(mp,pFlow0[k in brList;brData[k].zeroImpe == false], p0[k] == brData[k].g/(brData[k].τprime^2)*v0[brData[k].From]^2
               - brData[k].g/brData[k].τ*v0[brData[k].From]*v0[brData[k].To]*cos(θ0[brData[k].From] - θ0[brData[k].To] + brData[k].θtr)
               - brData[k].b/brData[k].τ*v0[brData[k].From]*v0[brData[k].To]*sin(θ0[brData[k].From] - θ0[brData[k].To] + brData[k].θtr)
               );
  @NLconstraint(mp,qFlow0[k in brList;brData[k].zeroImpe == false], q0[k] == (-brData[k].b - brData[k].bc/2)/(brData[k].τprime^2)*v0[brData[k].From]^2
               + brData[k].b/brData[k].τ*v0[brData[k].From]*v0[brData[k].To]*cos(θ0[brData[k].From] - θ0[brData[k].To] + brData[k].θtr)
               - brData[k].g/brData[k].τ*v0[brData[k].From]*v0[brData[k].To]*sin(θ0[brData[k].From] - θ0[brData[k].To] + brData[k].θtr)
               );
  @constraint(mp,zeroP0[k in brList;brData[k].zeroImpe == true], p0[k] == -p0[brData[k].revID]);
  @NLconstraint(mp,zeroQ0[k in brList;brData[k].zeroImpe == true], q0[k] + q0[brData[k].revID] == -brData[k].bc*v0[brData[k].To]^2);
  @constraint(mp,zeroV0[k in brList;brData[k].zeroImpe == true], v0[brData[k].To] == v0[brData[k].From]/brData[k].τ);
  @constraint(mp,zeroθ0[k in brList;brData[k].zeroImpe == true], θ0[brData[k].To] == θ0[brData[k].From] - brData[k].θtr);
  @constraint(mp,pBalance0[i in bList], sum(sp0[l] for l in bData[i].gen) == psh0[i] + bData[i].Pd + sum(p0[k] for k in brList if (brData[k].From == i)));
  @constraint(mp,qBalance0[i in bList], sum(sq0[l] for l in bData[i].gen) == qsh0[i] + bData[i].Qd + sum(q0[k] for k in brList if (brData[k].From == i)));

  # create the variables for the contingency cases
  @variable(mp,bData[i].Vmin <= v[i in bList, s in S] <= bData[i].Vmax);
  @variable(mp,gData[l].Pmin <= sp[l in gList, s in S;(!(l in contData[s].Loc))] <= gData[l].Pmax);
  @variable(mp,gData[l].Qmin <= sq[l in gList, s in S;(!(l in contData[s].Loc))] <= gData[l].Qmax);
  @variable(mp,p[k in brList, s in S]);
  @variable(mp,q[k in brList, s in S]);
  @variable(mp,psh[i in bList, s in S]);
  @variable(mp,qsh[i in bList, s in S]);
  @variable(mp,θ[i in bList, s in S]);
  @variable(mp,pΔ[s in S]);

  # create the constraints for the base case
  @NLconstraint(mp,flowBoundS[k in brList, s in S],p[k,s]^2 + q[k,s]^2 <= brData[k].t);
  @NLconstraint(mp,pShuntS[i in bList, s in S],psh[i,s] == bData[i].gsh*v[i,s]^2);
  @NLconstraint(mp,qShuntS[i in bList, s in S],qsh[i,s] == -bData[i].bsh*v[i,s]^2);
  @NLconstraint(mp,pFlowS[k in brList, s in S;(brData[k].zeroImpe == false)&(!(k in contDList[s].Loc))], p[k,s] == brData[k].g/(brData[k].τprime^2)*v[brData[k].From,s]^2
               - brData[k].g/brData[k].τ*v[brData[k].From,s]*v[brData[k].To,s]*cos(θ[brData[k].From,s] - θ[brData[k].To,s] + brData[k].θtr)
               - brData[k].b/brData[k].τ*v[brData[k].From,s]*v[brData[k].To,s]*sin(θ[brData[k].From,s] - θ[brData[k].To,s] + brData[k].θtr)
               );
  @NLconstraint(mp,qFlowS[k in brList, s in S;(brData[k].zeroImpe == false)&(!(k in contDList[s].Loc))], q[k,s] == (-brData[k].b - brData[k].bc/2)/(brData[k].τprime^2)*v[brData[k].From,s]^2
               + brData[k].b/brData[k].τ*v[brData[k].From,s]*v[brData[k].To,s]*cos(θ[brData[k].From,s] - θ[brData[k].To,s] + brData[k].θtr)
               - brData[k].g/brData[k].τ*v[brData[k].From,s]*v[brData[k].To,s]*sin(θ[brData[k].From,s] - θ[brData[k].To,s] + brData[k].θtr)
               );
  @constraint(mp,zeroPS[k in brList, s in S;brData[k].zeroImpe == true], p[k,s] == -p[brData[k].revID,s]);
  @NLconstraint(mp,zeroQS[k in brList, s in S;brData[k].zeroImpe == true], q[k,s] + q[brData[k].revID,s] == -brData[k].bc*v[brData[k].To,s]^2);
  @constraint(mp,zeroVS[k in brList, s in S;brData[k].zeroImpe == true], v[brData[k].To,s] == v[brData[k].From,s]/brData[k].τ);
  @constraint(mp,zeroθS[k in brList, s in S;brData[k].zeroImpe == true], θ[brData[k].To,s] == θ[brData[k].From,s] - brData[k].θtr);
  @constraint(mp,pBalanceS[i in bList, s in S], sum(sp[l,s] for l in bData[i].gen) == psh[i,s] + bData[i].Pd + sum(p[k,s] for k in brList if (brData[k].From == i)));
  @constraint(mp,qBalanceS[i in bList, s in S], sum(sq[l,s] for l in bData[i].gen) == qsh[i,s] + bData[i].Qd + sum(q[k,s] for k in brList if (brData[k].From == i)));
  @constraint(mp,recourse[l in gList, s in S], sp[l,s] == sp0[l] + gData[l].α*pΔ[s]);

  # create the Complementarity constraint
  @variable(mp, δplus[i in bList, s in S] >= 0);
  @variable(mp, δminus[i in bList, s in S] >= 0);
  @variable(mp, sqplus[i in bList,s in S] >= 0);
  @variable(mp, sqminus[i in bList,s in S] >= 0);
  @constraint(mp,vConstr1[i in bList,s in S], v0[i] - v[i,s] <= δplus[i,s]);
  @constraint(mp,vConstr2[i in bList,s in S], v[i,s] - v0[i] <= δminus[i,s]);
  @constraint(mp, spplusConstr[i in bList,s in S], sqplus[i,s] == sum(gData[l].Qmax - sq[l,s] for l in bData[i].gen));
  @constraint(mp, spminusConstr[i in bList,s in S], sqminus[i,s] == sum(sq[l,s] - gData[l].Qmin for l in bData[i].gen));
  for s in S
    for i in bList
      @complements(mp,0 <= δplus[i,s], sqplus[i,s] >= 0);
      @complements(mp,0 <= δminus[i,s], sqminus[i,s] >= 0);
    end
  end

  # add the contingency effect as constraints
  for s in S
    # for each contingency
    # first identify the type of the contingency
    if contData[s].Type == "B"
      # if it is a branch contingency
      for k in contData[s].Loc
        @constraint(mp, p[k,s] == 0.0);       # the real power flow through the branch is 0
        @constraint(mp, q[k,s] == 0.0);       # the reactive power flow through the branch is 0
      end
    elseif contData[s].Type == "T"
      # if it is a branch contingency
      for k in contData[s].Loc
        @constraint(mp, p[k,s] == 0.0);       # the real power flow through the branch is 0
        @constraint(mp, q[k,s] == 0.0);       # the reactive power flow through the branch is 0
      end
    elseif contData[s].Type == "G"
      # if it is a generator contingency
      for l in contData[s].Loc
        @constraint(mp, sp[l,s] == 0.0);
        @constraint(mp, sq[l,s] == 0.0);
      end
    end
  end

  # build the objective function
  @NLobjective(mp, Min, sum(sum(gData[l].cParams[n]*(sp0[l]*fData.baseMVA)^n for n in gData[l].cn) for l in gList));
  return mp;
end
