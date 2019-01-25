function d = dist_to_line(pt, v1, v2)
a = v1 - v2;
b = pt - v2;
d = norm(cross(a,b)) / norm(a);
if pt(2)<pt(1)
    d = -d;
end
end