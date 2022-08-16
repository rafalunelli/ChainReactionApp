function [image,sx,sy] = move_obj(image,sx,sy,vx,vy,dt)
    image.XData = image.XData + vx*dt;
    image.YData = image.YData + vy*dt;
    sx = sx + vx*dt;
    sy = sy + vy*dt;
end

