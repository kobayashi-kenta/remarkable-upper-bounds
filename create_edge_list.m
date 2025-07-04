%%  "edge_list(p+1, q+1)" denotes the index of edge        %%
%%  whoes mid point is  p/(2*n)*[1 0] + q/(2*n)*[a b] .    %%
function [edge_list, end_index] = create_edge_list(n, start_index)
    edge_list = zeros(2*n);
    i = start_index;
    for q = 0:n-1
        for p = 0:n-1-q
            cp = [p p+1 p];
            cq = [q q q+1];
            edge_list(cp(2)+cp(3)+1, cq(2)+cq(3)+1) = i;
            i = i + 1;
            edge_list(cp(3)+cp(1)+1, cq(3)+cq(1)+1) = i;
            i = i + 1;
            edge_list(cp(1)+cp(2)+1, cq(1)+cq(2)+1) = i;
            i = i + 1;
        end
    end
    end_index = i - 1;
end
