function [w] = clayton_cop1406(data, np, theta)

% Gera amostras da copula de Clayton
%RandStream.setDefaultStream(RandStream('mt19937ar','seed', sum(100*clock)));
w= zeros(np,size(data,2));
w(:,1)= rand(np,1);

for j=2:size(data,2)
  for i=1:np 
    %RandStream.setDefaultStream(RandStream('mt19937ar','seed', sum(100*clock)));
    t=rand;
    w(i,j)= ((sum(w(i,1:j-1).^-theta)-j+2).*(t.^(theta/(theta*(1-j)-1))-1) +1).^(-1/theta);
  end
end
end