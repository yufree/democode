t = [10, 30, 50, 70, 100]
r = [1, 5, 10, 25, 30, 50]
m = [1:10]
re = []
tindex = []
rindex = []
iindex = []

for i = 1:length(m)
   for k = 1:length(r)
        for j = 1:length(t)
            result = []
            for (l = 1:m(i))
                result = [result,1 / (m(l) ^ 2) * exp((-1020 * 9 * pi ^ 2 * m(l)^2 * t(j))/ (16 * r(k) ^ 2))]
                
            end
            re = [re,1 - 64 / (9 * pi ^ 2) * sum(result)]
            tindex = [tindex,t(j)]
            rindex = [rindex,r(k)]
            iindex = [iindex,m(i)]
        end
    end
end

d = [re;tindex;rindex;iindex]