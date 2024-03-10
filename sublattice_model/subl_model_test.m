% пример из системы Y-O, совпадает с ручными расчетами
s = sublattice_model('Y', 'O');
s = s.addSublattice(2, 'Y', 2, 'Y', 3).addSublattice(2, 'O', -2).addSublattice(1, 'O', -2, '', 0);
fprintf('\n\n%s\n\n', s.printFormula());
syms x_O
[comps, ys] = s.calculateModel('O', x_O);
for i = 1:size(comps, 1)
    fprintf('y(%s) = %s\n', comps(i), ys(i));
end

% пример из системы Fe-Cr
s = sublattice_model('Fe', 'Cr');
s = s.addSublattice(4, 'Cr', 0).addSublattice(16, 'Fe', 0, 'Cr', 0).addSublattice(10, 'Fe', 0);
fprintf('\n\n%s\n\n', s.printFormula());
syms x_Cr
[comps, ys] = s.calculateModel('Cr', x_Cr);
for i = 1:size(comps, 1)
    fprintf('y(%s) = %s\n', comps(i), ys(i));
end