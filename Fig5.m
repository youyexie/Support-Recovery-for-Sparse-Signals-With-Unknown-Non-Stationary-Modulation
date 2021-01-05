%% Exp-2

clear all;close all;clc


J = 3;
K = 3;
M = 150;
N = 100;
numIter = 100;
n_std =0.1;

%lambda = 2;
NsR = 0.02;
x = 3;
xlimit = 6;

PT = zeros(numIter,length(1:6));


    for J = 1:6
        %fprintf('eta = %d\n',eta);
        for Iter = 1:numIter
            fprintf('J= %d, Iter = %d',J,Iter);
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

            % Calculate the 2,infty norm of the recovery error
            deltaX_squared = max( norms( X-X0_scaled,2,1) )^2;
            fprintf(', error_squared = %f',deltaX_squared);
            
            if (length( find(norms( X, 2, 1 )>1e-5 )) == length( find((norms(X0_scaled,2,1)>1e-5)))  ) && min(  find(norms( X, 2, 1 )>1e-5 )==find((norms(X0_scaled,2,1)>1e-5)) )
                fprintf(', support recovery success.\n');
                PT(Iter,J) = deltaX_squared;
            else
                fprintf(', support recovery fail.\n');
            end
            
        end
    end


%%
figure(2);
%stdPT = std(PT);
PT(isnan(PT)) = 0;
meanPT = sum(PT)./sum(PT>0);

PT(PT==0)=NaN;
stdPT = nanstd(PT,[],1); 

plot(meanPT,'-','LineWidth',2);
hold on


%boxplot(PT)

for k2 = 1: length( meanPT )
   line([k2,k2],[meanPT(k2)+stdPT(k2),meanPT(k2)-stdPT(k2)],'LineWidth',2);%[max(logPT(:,k2)),min(logPT(:,k2))]);
   line([k2-0.1,k2+0.1],[meanPT(k2)+stdPT(k2),meanPT(k2)+stdPT(k2)],'LineWidth',2);%[max(logPT(:,k2)),max(logPT(:,k2))])
   line([k2-0.1,k2+0.1],[meanPT(k2)-stdPT(k2),meanPT(k2)-stdPT(k2)],'LineWidth',2);%[min(logPT(:,k2)),min(logPT(:,k2))])
end
scatter(1:6,meanPT,200,'+','r','LineWidth',3);

%title('Recovery error squared')
%yticks([1  3  5  7 9 11 13 15 ])
%yticklabels({'30','40','50','60','70','80','90','100'})
%xticks([1  2 3 4 5])
%xticklabels({'2','3','4','5','6'})
xlabel('J from 1 to 6') 
ylabel('Squared recovery error')
grid on
set(gca,'FontSize',20)
hold off
%%
%save('RecoveryError-J.mat','PT');
