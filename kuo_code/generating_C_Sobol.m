s = 6;
m = 4;

% Polynomial degree, coefficients a, and direction vector m
C(2).d = 1; C(2).a = [];  C(2).m = [1];
C(3).d = 2; C(3).a = 1; C(3).m = [1 1];
C(4).d = 3; C(4).a = 1; C(4).m = [1 3 7];
C(5).d = 3; C(5).a = 2; C(5).m = [1 1 5];
C(6).d = 4; C(6).a = 1; C(6).m = [1 3 1 1];



C(1).C = eye(m);
C(2).C = C(1).C;
for i = 2:m
    C(2).m = [C(2).m bitxor(2*C(2).m(end),C(2).m(end))];
    aa = (dec2bin(C(2).m(i))- '0')';
    empty = i - size(aa,1);
    C(2).C(1:i,i) = [zeros(empty,1); aa];
end

% Assuming a = a_d-1 * 2^0 + a_d-2 * 2^1 + ... + a_1 2^d-2
for j = 3:s
    C(j).C = C(1).C;
    for i = 2:m
        if i > size(C(j).m,2)
            vec_a = dec2bin(C(j).a)- '0';
            if size(vec_a,2) < C(j).d - 1
                vec_a = [zeros(1,C(j).d-size(vec_a,2)-1), vec_a];
            end
            vec_a = [(vec_a), 1, 1].*[2.^(1:C(j).d), 1]; % Matlab mistake: [fliplr(vec_a), 1, 1].*[2.^(1:C(j).d), 1]
            vec_a = vec_a.*[C(j).m(size(C(j).m,2):-1:size(C(j).m,2)-C(j).d+1), C(j).m(size(C(j).m,2)-C(j).d+1)];
            C(j).m = [C(j).m, bitxor_toni(vec_a)];
        end
        aa = (dec2bin(C(j).m(i))- '0')';
        empty = i - size(aa,1);
        C(j).C(1:i,i) = [zeros(empty,1);aa];
    end
end

x = zeros(2^m,s);
for n = 1:2^m-1
    vec_n = (dec2bin(n)- '0');
    vec_n = fliplr(vec_n)';
    for j = 1:s
        x(n+1,j) = 2.^-(1:m)*mod(C(j).C*[vec_n;zeros(m-size(vec_n,1),1)],2);
    end
end

sob = sobolset(s);
x_check = net(sob,2^m);
disp(x-x_check)

 