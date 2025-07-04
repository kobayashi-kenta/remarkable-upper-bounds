%%  "vertex_list(p+1, q+1)" denotes the index of vertex  %%
%%  whoes location is  p/n*[1 0] + q/n*[a b] .           %%
function [vertex_list, end_index] ...
        = create_vertex_list(n, start_index)
    vertex_list = zeros(n+1);
    i = start_index;
    for q = 0:n
        for p = 0:n-q
            vertex_list(p+1, q+1) = i;
            i = i + 1;
        end
    end
    end_index = i - 1;
end
