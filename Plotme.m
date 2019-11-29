function Plotme(t, time, x, xname, yname, zname, xlabelname, ylabelname, titlename, png_name)
    fig = figure();
    set(gcf,'color','w');
    x_t = subs(x, {t}, {time});
    plot(time, x_t(1:3, :), 'LineWidth', 2);
    grid on
    xlabel(xlabelname)
    ylabel(ylabelname)
    title(titlename)
    legend(xname, yname, zname, 'Fontsize', 14)
    frame = getframe(fig);
    im = frame2im(frame);
    [img,map] = rgb2ind(im,256);
    imwrite(img,map,png_name,'png');

end

