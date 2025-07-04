%%  "face_index(p+1, q+1)" denotes the index of face             %%
%%  whoes center of gravity is  p/(3*n)*[1 0] + q/(3*n)*[a b] .  %%
function [face_list, end_index] = create_face_list(n, start_index)
    face_list= zeros(3*n);
    i = start_index;
    for q = 0:n-1
        for p = 0:n-1-q
            cp = [p p+1 p];
            cq = [q q q+1];
            face_list(sum(cp)+1, sum(cq)+1) = i;
            i = i + 1;
            if p+q < n-1
	            cp = [p+1 p p+1];
    	        cq = [q+1 q+1 q];
                face_list(sum(cp)+1, sum(cq)+1) = i;
                i = i + 1;
            end
        end
    end
    end_index = i - 1;
end
