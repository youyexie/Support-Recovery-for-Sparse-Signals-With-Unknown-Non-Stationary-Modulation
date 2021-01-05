%% Exp-2

clear all;close all;clc


J = 3;
K = 3;
M = 150;
N = 100;
numIter = 50;
n_std =0.1;

%lambda = 2;
NsR = 0.02;
x = 3;

PT = zeros(15,20);
PT_NoZeroSolution = zeros(15,20);
PT_ExactRecovery = zeros(15,20);

for J = 1:20
    for N = 30:5:100
        %fprintf('eta = %d\n',eta);
        for Iter = 1:numIter
            fprintf('N = %d, J = %d, Iter = %d',N,J,Iter);
            % Dictionary A
            
            % Gaussian
            A = randn(N,M);
            
            % Subspace matrix B
            B = dftmtx(N)/sqrt(N);
            B = B(:,1:K);
            % Linear sensing process Phi
            Phi = zeros(N,K*M);
            for j=1:M
                for j2=1:K
                    Phi(:,K*(j-1)+j2) = diag(B(:,j2))*A(:,j);
                end
            end
            % Ground truth
            X0 = zeros(K,M);
            T = sort( randperm(M,J) );
            
            for k = 1:length(T)
                % h is unit norm
                %h = randn(K,1)+randn(K,1)*1i;
                %h = h/norm(h);
                X0(:,T(k)) = (randn(1)+randn(1)*1i)*(randn(K,1)+randn(K,1)*1i);
                
            end
                        
            % lambda
            gamma0 = sqrt( n_std^2*K*(log(M-J)+log(N)) );
            lambda = x*gamma0;
            
            % Observed signal
            NormX0 = norms(X0,2,1);
            minVecNorm = min( NormX0(NormX0>1e-5) );
            X0_scaled = (gamma0/NsR).*(X0./minVecNorm);
                                   
            NormX0_scaled = norms(X0_scaled,2,1);
            minVecNormX0_scaled = min( NormX0_scaled(NormX0_scaled>1e-5) );
            
            % Observed signal
            y = Phi*reshape(X0_scaled,[K*M,1]);
            
            % Noise
            n= n_std*randn([N,1])+n_std*randn([N,1])*1i;
            
            % Add noise to the observed signal
            y = y + n;
            

            
            cvx_begin quiet
            variable X(K,M) complex
            minimize( 0.5*pow_pos( norms( Phi*reshape(X,[K*M,1]) - y,2 ),2) + lambda*sum( norms( X, 2, 1 ) ) )
            cvx_end

            
            if (~ismember(0, ismember(find(norms( X, 2, 1 )>1e-5 ) , find((norms(X0_scaled,2,1)>1e-5)) ) ) )
                fprintf(', gamma0/minVecNorm :%f, support recovery success \n',gamma0/minVecNormX0_scaled)
                PT((N-30)/5+1,J) = PT((N-30)/5+1,J) + 1;
                % Only count the non-zero solution
                if (sum(norms(X))>1e-5)
                    PT_NoZeroSolution((N-30)/5+1,J) = PT_NoZeroSolution((N-30)/5+1,J) + 1;
                end
                 
                % Only count the exact support recovery
                if (length( find(norms( X, 2, 1 )>1e-5 )) == length( find((norms(X0_scaled,2,1)>1e-5)))  ) && min(  find(norms( X, 2, 1 )>1e-5 )==find((norms(X0_scaled,2,1)>1e-5)) )
                    PT_ExactRecovery((N-30)/5+1,J) = PT_ExactRecovery((N-30)/5+1,J) + 1;
                end
                
            else
                fprintf(', gamma0/minVecNorm :%f, support recovery fails\n',gamma0/minVecNormX0_scaled )
            end
            % Count the number of success recovery
            
        end
    end
    
end

%%
% figure(2);
% imagesc(PT/numIter);
% h=colorbar;
% title('0 Count')
% ylabel(h, 'Recovery rate')
% set(gca,'Ydir','normal')
% yticks([1  3  5  7 9 11 13 15 ])
% yticklabels({'30','40','50','60','70','80','90','100'})
% %xticks([1  3 5 7 9 11 13 15 17 19 21])
% %xticklabels({'5','7','9','11','13','15','17','19','21','23','25'})
% xlabel('J from 1 to 20') 
% ylabel('N from 30 to 100')
% 
% figure(3);
% imagesc(PT_NoZeroSolution/numIter);
% h=colorbar;
% title('0 Not Count')
% ylabel(h, 'Recovery rate')
% set(gca,'Ydir','normal')
% yticks([1  3  5  7 9 11 13 15 ])
% yticklabels({'30','40','50','60','70','80','90','100'})
% %xticks([1  3 5 7 9 11 13 15 17 19 21])
% %xticklabels({'5','7','9','11','13','15','17','19','21','23','25'})
% xlabel('J from 1 to 20') 
% ylabel('N from 30 to 100') 

figure(4);
imagesc(PT_ExactRecovery/numIter);
h=colorbar;
%title('0 Not Count')
ylabel(h, 'Exact support recovery rate')
set(gca,'Ydir','normal')
yticks([1  3  5  7 9 11 13 15 ])
yticklabels({'30','40','50','60','70','80','90','100'})
%xticks([1  3 5 7 9 11 13 15 17 19 21])
%xticklabels({'5','7','9','11','13','15','17','19','21','23','25'})
xlabel('J from 1 to 20') 
ylabel('N from 30 to 100') 
set(gca,'FontSize',20)
%%
%save('J.mat','PT');
%save('J(0NotCount).mat','PT_NoZeroSolution');
%save('J(ExactRecovery).mat','PT_ExactRecovery');