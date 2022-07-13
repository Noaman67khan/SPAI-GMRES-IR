function M = spai_ss(A,espai,alpha,beta)
% SPAI_SS Constructs a sparse approximate inverse in single precision
%   A is the input matrix
%   espai is the tolerance parameter for the quality of each column
%   alpha is the number of steps for each column
%   beta is the number of steps for adding nonzeros to the pattern
%   M is the output sparse approximate inverse

A = single(full(A));
n = length(A);
J = spones(A');
I = eye(n);
M = single(zeros(n));
for k = 1:n
    %ek is kth column of identity
    ek = I(:,k);
    %Get index set Jk
    Jk = find(J(:,k));
  
    for step = 1:alpha
        %Construct set Ik
        Ik=[];
        for i=1:n
            if sum(abs(A(i,Jk)))~=0
                Ik=[Ik i];
            end
        end
        
        Atk = single(A(Ik,Jk));
        etk = ek(Ik);
        [Qt,Rt] = qr(Atk,0);

        Mtk = single(Rt)\(single(Qt')*single(etk));
        rtk = single(Atk)*single(Mtk) - single(etk);
        if norm(rtk) < espai
            %[norm(rtk),k];
            break
        end
         
       
        %Construct index set Lk
        Lk = union(Ik,k);
        
        Jtk=[];
        for ll = 1:numel(Lk)
            l=Lk(ll);
            Nl=[];
            for j=1:n
                if A(l,j)~=0
                    Nl = union(Nl,j);
                end
            end
            Jtk = union(Jtk, Nl);
        end
        Jtk = setdiff(Jtk,Jk);
        
        rok = 0;
        Rojk = [];
        n1 = norm(rtk);
       
        for jj = 1:numel(Jtk)
            j = Jtk(jj);
            
            n2 = norm(full(A(Ik,j)));
          
            rojk = abs((n1^2-((rtk'*A(Ik,j))^2/(n2^2))))^0.5;
            rok = rok + rojk;
            newrow = [j rojk];
            Rojk = [Rojk ; newrow];
        end
        
        rok = rok/(numel(Jtk));
        
        for idx=1:beta
            [min_val,min_idx]= min(Rojk(:,2));
            j = Rojk(min_idx,1);
            if(Rojk(min_idx,2) <= rok)
                Jk = union(Jk,j);
                Jtk = setdiff(Jtk,j);
                Rojk(min_idx,:) = [];
            end
        end
        
        
    end
    %disp(k)
    M(Jk,k)= Mtk;    
end
end