classdef dirichlet < handle
    properties
        dim
        alpha_0
        alpha
    end
    
    methods
        function self = dirichlet(dim,alpha_0,alpha)
            self.dim = dim;
            
            if ~exist('alpha_0','var')
                self.alpha_0 = ones(self.dim,1);
            else
                self.alpha_0 = alpha_0;
            end
            
            if ~exist('alpha','var')
                self.alpha = self.alpha_0.*(1+rand(size(self.alpha_0)));
            else
                self.alpha = alpha;
            end
        end
        
        function update(self,data,beta)
            if ~exist('data','var')
                data = 0;
            end
            if ~exist('beta','var')
                beta = 1;
            end
            self.alpha = beta .* self.alpha_0 + data;
            if(sum(isnan(self.alpha))>0)
                'NaNs detected'
                self.alpha=self.alpha_0;
            end
        end
        
        function updateSS(self,NA)
            if ~exist('NA','var')
                NA = 0;
            end
            
            self.alpha = self.alpha_0 + NA;
            if(sum(isnan(self.alpha))>0)
                'NaNs detected'
                self.alpha=self.alpha_0;
            end
        end
                
        function rawupdate(self,data,p)
            if(~exist('p','var'))
                p=ones(1,size(data,2));
            end
            idx=find(~isnan(sum(data,2)));
            SEx = data(:,idx)*p(idx)';
            self.update(SEx); 
        end
        
        function res = mean(self)
            res = self.alpha / sum(self.alpha(:));
        end
        
        function res = geomean(self)
            res = exp(self.loggeomean());
        end
        
        function res = loggeomean(self)
            res = bsxfun(@minus,psi(self.alpha),psi(sum(self.alpha(:))));
        end

        function res = variance(self)
           alpha_sum = sum(self.alpha(:));
           res = self.alpha.*(alpha_sum-self.alpha)/alpha_sum^2/(alpha_sum-1);
        end
        
        function res = covariance(self)
            alpha_sum = sum(self.alpha(:));
            res = -self.alpha'*self.alpha/alpha_sum^2/(alpha_sum-1);
            res = res.*(ones(self.dim)-eye(self.dim))+diag(self.variance);
        end
        
        function res = entropy(self)
            alpha_sum = sum(self.alpha(:));
            res = sum(gammaln(self.alpha(:))) - gammaln(alpha_sum) + ...
                (alpha_sum - self.dim)*psi(alpha_sum) - ...
                sum((self.alpha(:) - 1).*psi(self.alpha(:)));
        end

        function res = KLqprior(self)
            alpha_sum = sum(self.alpha);
            res = gammaln(alpha_sum) - sum(gammaln(self.alpha)) ...
                - gammaln(sum(self.alpha_0)) + sum(gammaln(self.alpha_0)) ...
                + sum((self.alpha-self.alpha_0).*(psi(self.alpha)-psi(alpha_sum)));
        end

        function res = Eloglikelihood(self,data) % assumes multinomial observations
            res = data*self.loggeomean + gammaln(1+sum(data,2)) - sum(gammaln(1+data),2);
            res(isnan(res)) = 0;
        end
        
        function res = expectlogjoint(self,beta)
            if ~exist('beta','var')
                alpha_prior = self.alpha_0;
            else
                alpha_prior = self.alpha_0 .* beta;
            end
            
            contrib = (alpha_prior - 1).*self.loggeomean();
            
            res =  - sum(gammaln(alpha_prior(:))) + ...
                gammaln(sum(alpha_prior(:))) + ...
                sum(contrib(:));
        end
        
        function res = lowerboundcontrib(self,beta)
            if ~exist('beta','var')
               beta = 1;
            end
            res = self.entropy() + self.expectlogjoint(beta);
        end
    end
end
