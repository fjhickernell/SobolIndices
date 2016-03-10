m = 6;

aux = eye(m);
for i = 0:2^(m*(m-1)/2)-1
    a = dec2bin(i) - '0';
    empty = m*(m-1)/2 - size(a,2);
    if empty
        a = [zeros(1,empty) a];
    end
    ind = [];
    for j = 2:m
        ind = [ind j+m*(j-2):m*(j-1)];
    end
    A = aux;
    A(ind) = a;
    A = A';
    for k = 0:m-2
        if mod(A^(2^k),2) == aux
            GROUP.order(i+1) = 2^k;
            GROUP.element(i+1).matrix = A;
            break
        elseif k == m-2
            GROUP.order(i+1) = 2^(m-1);
            GROUP.element(i+1).matrix = A;
        end
    end
end

disp(max(GROUP.order))
    
    