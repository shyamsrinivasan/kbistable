function [time,dGene,ngenes,ntf,model] = testlinear(P,TF,DNA,Gene1,PGenes,EGene,EPSGene,M,model)
%linlogfunc1 - Highest accuracy/approximation closest to reality
      
    ngenes = size(Gene1,1);
    ntf = size(TF,1);
    
    KB = P.basaltr;         %Basal Expression
    KS = P.regtr;           %Regulated Expression
    beta = P.actcoeff;             
    theta = P.repcoeff;            
    n = P.hillcoeff;                
    beta1 = P.actcoeff;            %Metabolite Activation
    theta1 = P.repcoeff;          %Metabolite Repression
    
    ktranslation = P.ktranslation;    
    
    KD = P.genedecay;
    protdecay = P.protdecay;
    
    V = [Gene1;TF];
    Vinit = V';
   
    flag=zeros(ngenes,1);
    for k=1:ngenes
        if ~isempty(find(EGene(k,:)>0,1))
            flag(k) =1;
        end
        if ~isempty(find(EGene(k,:)<0,1))
            flag(k)=2;
        end
        if ~isempty(find(EGene(k,:)<0,1)) && ~isempty(find(EGene(k,:)>0,1))
            flag(k)=3;
        end
        if ~isempty(find(EPSGene(k,:)>0,1))
            flag(k) =4;
        end
        if ~isempty(find(EPSGene(k,:)<0,1))
            flag(k)=5;
        end
        if ~isempty(find(EPSGene(k,:)<0,1)) && ~isempty(find(EPSGene(k,:)>0,1))
            flag(k)=6;
        end
        if ~isempty(find(EGene(k,:)>0,1)) && ~isempty(find(EPSGene(k,:)>0,1))
            flag(k)=7;
        end
        if ~isempty(find(EGene(k,:)>0,1)) && ~isempty(find(EPSGene(k,:)<0,1))
            flag(k)=8;
        end
        if ~isempty(find(EGene(k,:)<0,1)) && ~isempty(find(EPSGene(k,:)>0,1))
            flag(k)=9;
        end
        if ~isempty(find(EGene(k,:)<0,1)) && ~isempty(find(EPSGene(k,:)<0,1))
            flag(k)=10;
        end
        if ~isempty(find(EGene(k,:)<0,1)) && ~isempty(find(EGene(k,:)>0,1)) && ~isempty(find(EPSGene(k,:)<0,1))
            flag(k)=11;
        end
        if ~isempty(find(EGene(k,:)<0,1)) && ~isempty(find(EGene(k,:)>0,1)) && ~isempty(find(EPSGene(k,:)>0,1))
            flag(k)=12;
        end
        if ~isempty(find(EGene(k,:)>0,1)) && ~isempty(find(EPSGene(k,:)>0,1)) && ~isempty(find(EPSGene(k,:)<0,1))
            flag(k)=13;
        end
        if ~isempty(find(EGene(k,:)<0,1)) && ~isempty(find(EPSGene(k,:)>0,1)) && ~isempty(find(EPSGene(k,:)<0,1))
            flag(k)=14;
        end
    end

    model.flag = flag;
    [time, dG] = ode23(@linlogfunc1,[0 72000],Vinit);
    dGene = dG;
        
    function dV = linlogfunc1(t,V)
       dV = ones(ngenes+ntf,1);
       G = V(1:ngenes);
       TF = V(ngenes+1:end);
        
       for i=1:ngenes
           
          switch flag(i)
              case 0%No R
                  dV(i) = KB(i)*(DNA(i))-kd*V(i);
          
              case 1%ER>0
                  dV(i) = KS(i)*sum((TF(EGene(i,:)>0))./...
                          (beta(i)+TF(EGene(i,:)>0)))-KD*V(i);
                  
              case 2%ER<0
                  dV(i) = KS(i)*sum(1/...
                          (1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))-KD*V(i);
                  
              case 3%ER>0 && ER<0
                  dV(i) = KS(i)*(sum(1/(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))+...
                          sum((TF(EGene(i,:)>0)*beta(i))./(beta(i)+TF(EGene(i,:)>0))))...
                          -KD*V(i);
                  
              case 4%MR>0
                  dV(i) = KS(i)*sum((M(EPSGene(i,:)>0))./(beta1(i)+M(EPSGene(i,:)>0)))-...
                          KD*V(i);
                  
              case 5%MR<0
                  dV(i) = KS(i)*(sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i))))-...
                          KD*V(i);%M(EPSGene(i,:)<0)
                  
              case 6%MR>0 && MR<0
                  dV(i) = KS(i)*(sum((M(EPSGene(i,:)>0))./(beta1(i)+M(EPSGene(i,:)>0)))+...
                          sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i))))-...
                          KD*V(i);%M(EPSGene(i,:)<0).
                  
              case 7%ER>0 && MR>0
                  dV(i) = KS(i)*(sum((TF(EGene(i,:)>0))./(beta(i)+TF(EGene(i,:)>0))))*...
                          (sum((M(EPSGene(i,:)>0))./(beta1(i)+M(EPSGene(i,:)>0))))-...
                          KD*V(i); 
                  
              case 8%ER>0 && MR<0
                  dV(i) = KS(i)*(sum((TF(EGene(i,:)>0))./(beta(i)+TF(EGene(i,:)>0))))*...
                          (sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i))))-...
                          KD*V(i);%M(EPSGene(i,:)<0).
                  
              case 9%ER<0 && MR>0
                  dV(i) = KS(i)*sum(1./(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))*...
                          sum((M(EPSGene(i,:)>0))./(beta1(i)+M(EPSGene(i,:)>0)))-...
                          KD*V(i);
                  
              case 10%ER<0 && MR<0
                  dV(i) = KS(i)*sum(1/(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))*...
                          sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i)))-...
                          KD*V(i);%TF(EGene(i,:)<0). M(EPSGene(i,:)<0).
                  
              case 11%ER>0 && ER<0 && MR<0
                  dV(i) = KS(i)*(sum(1./(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))+...
                          sum((TF(EGene(i,:)>0))./(beta(i)+TF(EGene(i,:)>0))))*...
                          sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i)))-...
                          KD*V(i);%M(EPSGene(i,:)<0).
                  
              case 12%ER>0 && ER<0 && MR>0
                  dV(i) = KS(i)*(sum(1./(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))+...
                          sum((TF(EGene(i,:)>0))./(beta(i)+TF(EGene(i,:)>0))))*...
                          sum((M(EPSGene(i,:)>0))./(beta1(i)+M(EPSGene(i,:)>0)))-...
                          KD*V(i);
                  
              case 13%ER>0 && MR<0 && MR>0
                  dV(i) = KS(i)*(sum((M(EPSGene(i,:)>0))./(beta1(i)+M(EPSGene(i,:)>0)))+...
                          sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i))))*...
                          sum((TF(EGene(i,:)>0))./(beta(i)+TF(EGene(i,:)>0)))-...
                          KD*V(i);%M(EPSGene(i,:)<0).
                  
              case 14%ER<0 && MR<0 && MR>0
                  dV(i) = KS(i)*(sum((M(EPSGene(i,:)>0))./(beta1(i)+M(EPSGene(i,:)>0)))+...
                          sum(1/(1+(M(EPSGene(i,:)<0)/theta1(i)).^n(i))))*...
                          sum(1./(1+(TF(EGene(i,:)<0)/theta(i)).^n(i)))-...
                          KD*V(i);%M(EPSGene(i,:)<0).
          end
       end
       
       %Translation
       dV(ngenes+1:end) = ktranslation*V(1:ngenes)-0.0006*V(ngenes+1:end);
       
    end

    function dV = linlogfunc(t,V)
       dV = ones(ngenes+ntf,1);
       G = V(1:ngenes);
       TF = V(ngenes+ntf:end);
        
       for i=1:ngenes
                      
          Ractivated = 0;
          Rinhibited=0;
          if ~isempty(find(EGene(i,:)<0,1))
             Rinhibited = sum(1./(TF(EGene(i,:)<0)));
          end
          if ~isempty(find(EGene(i,:)>0,1))
             Ractivated = sum(TF(EGene(i,:)>0));
          end
           
          Ts =1;
          if any(EPSGene(i,:)>0)||any(EPSGene(i,:)<0)
              if any(EPSGene(i,:)<0)
                 Ts = Ts*~any(M(EPSGene(i,:)<0));
              else                              % => any(EPSGene(i,:)>0)
                 Ts = Ts*any(M(EPSGene(i,:)>0));
              end
          end
                       
          dV(i) = kb*(DNA(i))+ks*(Ractivated+Rinhibited)*Ts-protdecay*V(i);
              
                    
       end
       dV(ngenes+1:end) = ktranslation*G-protdecay*TF;
    end        
  
end        
