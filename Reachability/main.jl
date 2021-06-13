using LazySets

#Z = Zonotope([1.0, 0.0] , [0.1 0.0; 0.0 0.1])

A = [-1.0 -4.0; 4.0 -1.0]

#A_Exp = exp(A*0.01)

#using LinearMaps
#linear_map(A_exp, Z)

#vertices_list(Z)

f = open("zonotopes.txt", "r") 
out = open("zonotopes_vertices.txt","w")
  
while ! eof(f)   
    s = readline(f)           
    s = split(s, " ")

    n = parse(Int, s[1])
    c = Vector{Float64}(undef, n)
    for i = 2:n+1
        val = parse(Float64, s[i])
        c[i-1] = val
    end

    s = readline(f)           
    s = split(s, " ")

    m = parse(Int, s[1])
    g = Matrix{Float64}(undef, m, n)
    for i = 1:m
        for j = 1:n
            val = parse(Float64, s[1 + (i-1)*n + j])
            g[i, j] = val
        end
    end

    Z = Zonotope(c, g)

    v = vertices_list(Z)
    n = size(v)[1]

    write(out, string(size(v)[1]))
    write(out, " ")
    for i=1:n
        print(v[i])
        for j=1:2
            write(out, string(v[i][j]))
            write(out, " ")
        end
    end
    write(out, "\n")
 end 

close(f) 
close(out)
