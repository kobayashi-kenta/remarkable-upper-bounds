function [triangle, orientation] ...
        = set_triangle_beta(n, vertex_list, edge_list)
    triangle = zeros(n^2, 6);
    orientation = zeros(n^2, 1);
    i = 0;
    for q = 0:n-1
        for p = 0:n-1-q
            cp = [p p+1 p];
            cq = [q q q+1];
            vertex1  = vertex_list(cp(1)+1, cq(1)+1);
            vertex2  = vertex_list(cp(2)+1, cq(2)+1);
            vertex3  = vertex_list(cp(3)+1, cq(3)+1);
            edge1 = edge_list(cp(2)+cp(3)+1, cq(2)+cq(3)+1);
            edge2 = edge_list(cp(3)+cp(1)+1, cq(3)+cq(1)+1);
            edge3 = edge_list(cp(1)+cp(2)+1, cq(1)+cq(2)+1);
            i = i + 1;
            triangle(i, :) ...
            = [vertex1, vertex2, vertex3, edge1, edge2, edge3];
            orientation(i) = 1;
            if p+q < n-1
	            cp = [p+1 p p+1];
    	        cq = [q+1 q+1 q];
                vertex1  = vertex_list(cp(1)+1, cq(1)+1);
                vertex2  = vertex_list(cp(2)+1, cq(2)+1);
                vertex3  = vertex_list(cp(3)+1, cq(3)+1);
                edge1 = edge_list(cp(2)+cp(3)+1, cq(2)+cq(3)+1);
                edge2 = edge_list(cp(3)+cp(1)+1, cq(3)+cq(1)+1);
                edge3 = edge_list(cp(1)+cp(2)+1, cq(1)+cq(2)+1);
                i = i + 1;
                triangle(i, :) ...
                = [vertex1, vertex2, vertex3, edge1, edge2, edge3];
                orientation(i) = 2;
            end
        end
    end
end
