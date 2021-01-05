%% Exp-2

clear all;close all;clc


J = 3;
K = 3;
M = 150;
N = 100;
numIter = 50;
n_std =0.1; %0.5

%lambda = 2;%5

gamma0 = sqrt( n_std^2*K*(log(M-J)+log(N)) );

xlimit = 6;
PT = zeros(length(0.02:0.02:1),length(0.2:0.2:xlimit));
PT_NoZeroSolution = zeros(length(0.02:0.02:1),length(0.2:0.2:xlimit));
PT_ExactRecovery = zeros(length(0.02:0.02:1),length(0.2:0.2:xlimit));

for NsR = 0.02:0.02:1
    for x = 0.2:0.2:xlimit
        %fprintf('eta = %d\n',eta);
        for Iter = 1:numIter
            fprintf('NsR = %f, x = %f, Iter = %d',NsR,x,Iter);
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
                temp = (randn(1)+randn(1)*1i)*(randn(K,1)+randn(K,1)*1i);
                X0(:,T(k)) = temp; % Try not control the ratio
                
            end
            % Observed signal
            NormX0 = norms(X0,2,1);
            minVecNorm = min( NormX0(NormX0>1e-5) );
            X0_scaled = (gamma0/NsR).*(X0./minVecNorm);
            
%             % Testing
%             NormX0_scaled = norms(X0_scaled,2,1);
%             test_minVecNorm = min( NormX0_scaled(NormX0_scaled>1e-5) )
            
            y = Phi*reshape(X0_scaled,[K*M,1]);
            
            % Noise
            n= n_std*randn([N,1])+n_std*randn([N,1])*1i;
            
            % Add noise to the observed signal
            y = y + n;
            
            % Lambda
            lambda = x*gamma0;
            
            cvx_begin quiet
            variable X(K,M) complex
            minimize( 0.5*pow_pos( norms( Phi*reshape(X,[K*M,1]) - y,2 ),2) + lambda*sum( norms( X, 2, 1 ) ) )
            cvx_end

            % Count if the support of the solution is a subset of the
            % ground truth
            if (~ismember(0, ismember(find(norms( X, 2, 1 )>1e-5 ) , find((norms(X0_scaled,2,1)>1e-5)) ) ))
                fprintf(', support recovery success\n')
                PT(int32(NsR*50),int32(x*5)) = PT(int32(NsR*50),int32(x*5)) + 1;
               % Only count the non-zero solution
                if length( find(norms( X, 2, 1 )>1e-5 ))>=1
                    PT_NoZeroSolution(int32(NsR*50),int32(x*5)) = PT_NoZeroSolution(int32(NsR*50),int32(x*5)) + 1;
                end
                
                % Only count the exact support recovery
                if (length( find(norms( X, 2, 1 )>1e-5 )) == length( find((norms(X0_scaled,2,1)>1e-5)))  ) && min(  find(norms( X, 2, 1 )>1e-5 )==find((norms(X0_scaled,2,1)>1e-5)) )
                    PT_ExactRecovery(int32(NsR*50),int32(x*5)) = PT_ExactRecovery(int32(NsR*50),int32(x*5)) + 1;
                end
            else
                fprintf(', support recovery fails\n')
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
% %yticks([1  3  5  7 9 11 13 15 ])
% yticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
% %xticks([1  3 5 7 9 11 13 15 17 19 21])
% xticklabels({'1','2','3','4','5','6'})
% xlabel('x from 1 to 6') 
% ylabel('NsR from 0.1 to 1') 
% 
% figure(3)
% imagesc( PT_NoZeroSolution/numIter);
% h=colorbar;
% title('0 Not Count')
% ylabel(h, 'Recovery rate')
% set(gca,'Ydir','normal')
% %yticks([1  3  5  7 9 11 13 15 ])
% yticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
% %xticks([1  3 5 7 9 11 13 15 17 19 21])
% xticklabels({'1','2','3','4','5','6'})
% xlabel('x from 1 to 6') 
% ylabel('NsR from 0.1 to 1') 

figure(4)
imagesc( PT_ExactRecovery/numIter);
h=colorbar;
%title('Exact Recovery')
ylabel(h, 'Exact support recovery rate')
set(gca,'Ydir','normal')
%yticks([1  3  5  7 9 11 13 15 ])
yticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
%xticks([1  3 5 7 9 11 13 15 17 19 21])
xticklabels({'1','2','3','4','5','6'})
xlabel('k from 0.2 to 6') 
ylabel('\gamma  from 0.02 to 1') 
set(gca,'FontSize',15)
%%
%save('lambda','PT');
%save('lambda(0NotCount).mat','PT_NoZeroSolution');
%save('lambda(ExactRecovery).mat','PT_ExactRecovery');