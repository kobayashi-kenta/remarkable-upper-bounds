function triangle = set_triangle_alpha(n, edge_list, face_list)
    triangle = zeros(n^2, 4);
    i = 0;
    for q = 0:n-1
        for p = 0:n-1-q
            cp = [p p+1 p];
            cq = [q q q+1];
            edge1 = edge_list(cp(2)+cp(3)+1, cq(2)+cq(3)+1);
            edge2 = edge_list(cp(3)+cp(1)+1, cq(3)+cq(1)+1);
            edge3 = edge_list(cp(1)+cp(2)+1, cq(1)+cq(2)+1);
            face  = face_list(sum(cp)+1, sum(cq)+1);
            i = i + 1;
            triangle(i, :) = [edge1, edge2, edge3, face];
            if p+q < n-1
                cp = [p+1 p p+1];
                cq = [q+1 q+1 q];
                edge1 = edge_list(cp(2)+cp(3)+1, cq(2)+cq(3)+1);
                edge2 = edge_list(cp(3)+cp(1)+1, cq(3)+cq(1)+1);
                edge3 = edge_list(cp(1)+cp(2)+1, cq(1)+cq(2)+1);
                face  = face_list(sum(cp)+1, sum(cq)+1);
                i = i + 1;
                triangle(i, :) = [edge1, edge2, edge3, face];
            end
        end
    end
end
