function [school_pop, Phi] = SchoolPopulation(YearSize, YearGroups, NumCloseContacts, p_rewire)

%This function makes the school population matrix of year groups for
%the longer version of the code

%N.B. Currently all year groups are assumed to have the same close contact
%structure, but could be adjusted so that this was not the case

%N.B. p_rewire can be used to make non-exclusive close contact groups, this
%is not used in the paper (i.e. we set p_rewire = 0)- for p_rewire > 0, the
%model sometimes ends up in an infinite loop

if nargin == 0
    YearSize = 200;
    YearGroups = 5;
    NumCloseContacts = 25; % must be a divisor of the relevant yearsize
    p_rewire = 0;
end



NumGroups = YearSize/NumCloseContacts; %number of close contact groups in each year group

school_matrix = [];
yeargroup_vec = [];

for i = 1:YearGroups
    
%Make year group population    
    M = [];

for j = 1:NumGroups
    
    M = blkdiag(M, sparse(ones(NumCloseContacts)));
    
end

%remove self contacts
M = M - speye(length(M));
M = triu(M);
[x,y]=find(M);   % links without the diagonal

r = randperm(length(x), round(length(x)*p_rewire)); % randomly choose links to swap

%making the length of r even if r is odd
if mod(length(r),2) ~= 0
    r(end) = [];
end


    while ~isempty(r)
        
        X = x(r(1)); Y = y(r(1));
                
        if M(X,Y) == 1 %if the link is still there      
        
            nX = 1; nY = 1;
            % this while loop checks that all 4 people are unique
            %this can go on infinitely...
            while length(unique([X Y nX nY]))<4  || M(X,nY) + M(nY,X) + M(nX,Y) + M(Y,nX) >0
                r2 = randi(length(r) - 1) + 1;
               
                nX = x(r(r2)); nY = y(r(r2));     
            end
            
            M(X,Y) = 0; %remove link
            M(nX,nY) =0; %remove link
            
            if X < nY %so that we add link to upper diagonal matrix
                M(X,nY) = 1; %add link
            else
                M(nY, X) = 1;
            end
            
            if nX < Y
                M(nX, Y) = 1; %add link
            else
                M(Y, nX) = 1;
            end
            
            %deleting used r elements
            %need to delete this one first as the other way round would
            %change the element that r2 refers to
            r(r2) = [];
            r(1) = [];
            
        end
        
    end
    
    M = M + M';
   % Triples = sum(sum(M2)) - trace(M2);
   % Triangles = trace(M2*M);
   % Phi(i) = Triples/Triangles;
   Phi(i) = 0; %Phi not used so currently not calculates
    
    school_matrix = blkdiag(school_matrix, M); %store school matrix
    yeargroup_vec = [yeargroup_vec, i*ones(1,YearSize)]; %store which year group each individual is in
    
    
end

%store information in school_pop struct
school_pop.school_matrix = school_matrix;
school_pop.yeargroup_vec = yeargroup_vec;
school_pop.yearsize = YearSize;
school_pop.yeargroups = YearGroups;