function [school_pop] = SchoolPopulationquicker(YearSize, YearGroups, NumCloseContacts)

%This function makes the school population matrix of year groups for
%the quicker version of the code

%N.B. Only one school year group matrix is made - assumed all year groups
%have the same close contact structure

if nargin == 0
    YearSize = 200;
    YearGroups = 5;
    NumCloseContacts = 25; % must be a divisor of the relevant yearsize
    p_rewire = 0;
end


NumGroups = YearSize/NumCloseContacts; %number of close contact groups in a year

yeargroup_vec = []; %not stored in this version

%create school_matrix
A = sparse(ones(NumCloseContacts));
A = A - speye(length(A));
Ac = repmat({A}, 1, NumGroups);    
school_matrix = blkdiag(Ac{:});

%store school_matrix, as well as YearGroups and YearSize
school_pop.school_matrix = school_matrix;
school_pop.yeargroup_vec = yeargroup_vec;
school_pop.yeargroups = YearGroups;
school_pop.yearsize = YearSize;
