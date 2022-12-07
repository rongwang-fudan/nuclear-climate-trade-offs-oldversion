function x2 = timeconvert( realtime, x, yearstart, int )
%   Author:  Rong Wang // contact rongwang@fudan.edu.cn //
%   Date:  2021.7.12
%   Conversion of a time series to the time format of realtime

T = size(realtime,1);
[T2,n] = size(x);
x1 = ones(T,n).*(-999);

i=1;
j=max(1,floor((realtime(i,1)-yearstart)/int)+1);
while j<=T2 && i<=T
    x1(i,1:n)=x(j,1:n);
    i=i+1;
    j=max(1,floor((realtime(i,1)-yearstart)/int)+1);
end

x2 = x1(1:(i-1),:);

end

