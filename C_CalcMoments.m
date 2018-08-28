function moments = C_CalcMoments(Input,Ps,Es,varargin)
% Estimate the 1st and 2nd moments of the result of some test function

% Update online if necessary
if(nargin>3) [~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'BfFields',[1,2]);


% If we get state data, run the test on each state to form a history
if(size(Input,3)>1)     
    bfhist(1:size(Input,3),1) = 1:size(Input,3);
    for ii=1:size(Input,3)
        temp = Es.TestFunc(Input(:,:,ii),Ps,Es);
        bfhist(ii,2:1+length(temp))=temp;
    end;
else          % Or, assume we got a history
    bfhist = Input;
end;

if(size(bfhist,1)>1)
    dt = diff(bfhist(1:2,Es.BfFields(1)));

    moments(1) = dt*sum(1-bfhist(:,Es.BfFields(2)));
    moments(2) = dt*sum((1-bfhist(:,Es.BfFields(2))).^2);
else
    moments = [NaN NaN];
end;
end


