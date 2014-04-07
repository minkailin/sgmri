pro split, ntot=ntot, kmin=kmin, kmax=kmax, kmid=kmid

dk = kmax - kmin

imid = (kmid - kmin)*(ntot - 1d0)/dk + 1

imid = fix(imid)

kmid = kmin + dk*(imid -1d0)/(ntot-1d0)

kmidp1=kmin + dk*(imid)/(ntot-1d0)

print,'imid    =', imid
print, format='("kmid    =",e22.15)', kmid
print, format='("kmidp1  =",e22.15)', kmidp1
print,'n2      =', ntot - imid 

end
