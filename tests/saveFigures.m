function saveFigures(fName, pos)
    set(gcf, 'PaperPositionMode', 'auto');
    if nargin > 1
        set(gcf, 'pos', pos);
    end
    print(gcf,'-dpng', '-noui', '-loose', fName);
end