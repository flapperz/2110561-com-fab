// this will be overwrited by auto generate city
unit_d = 100;

function side_range(w) = [-((w-2)/2):(w-2)/2];

module window()
{   
    difference()
    {
        cube([40,10,60],center=true);
        cube([30,10+$fs,50],center=true);
    }
    
    cube([40,15,5],center=true);
    translate([0,0,-30])
        cube([40,15,5], center=true);
}

module tower(H)
{
    cylinder(r=30, h=100*H);
    cylinder(r=40, h=10);
    translate([0,0,100*H])
        difference()
        {
            cylinder(r=45, h=40);
            translate([0,0,20+$fs]) cylinder(r=35, h=20);
            translate([0,0,35])
            for (i=[0:5])
                rotate([0,0,60*i]) cube([100,20,10+$fs], center=true);
        }
}

module wall(H)
{
    difference()
    {
        translate([0,0,H/2*100])
            cube([100,20,H*100],center=true);
        
        translate([0,0,0])
            cube([30,100,70],center=true);
    }
    
    for (i=[0:1])
        translate([-i*50,-10,H*100])
            cube([25,20,20]);
}

module wall_edge(H)
{
    translate([25,0,H/2*100])
    cube([50,20,H*100],center=true);

    rotate([0,0,90])
        translate([25,0,H/2*100])
            cube([50,20,H*100],center=true);

    translate([-10,-10,0])
        cube([20,20,H*100]);
}

module house1()
{
    roof_r = 40;

    translate([0,0,40]) cube(80,center=true);
    translate([0,0,80])
        polyhedron(
            points=[ [roof_r,roof_r,0],[roof_r,-roof_r,0],
                    [-roof_r,-roof_r,0],[-roof_r,roof_r,0], // the four points at base
                    [0,0,roof_r]  ],                                 // the apex point 
            faces=[ [0,1,4],[1,2,4],[2,3,4],[3,0,4],              // each triangle side
                    [1,0,3],[2,1,3] ]                         // two triangles for square base
    );

    for (i=[0:3])
        rotate([0,0,i*90])
        translate([0,-40,40])
            window();   
}

module house2()
{   
    roof_r = 35;

    for (i=[0:3])
        rotate([0,0,i*90])
        translate([30,30,0])
            cylinder(r=10, h=70);

    translate([0,0,5])cube([90,90,10],center=true);
    translate([0,0,70]) cube([90,90,10],center=true);
    
}

module house3()
{
    roof_r = 40;

    translate([0,0,80]) cube([80,80,160],center=true);
    translate([0,0,160])
        polyhedron(
            points=[ [roof_r,roof_r,0],[roof_r,-roof_r,0],
                    [-roof_r,-roof_r,0],[-roof_r,roof_r,0], // the four points at base
                    [0,0,roof_r]  ],                                 // the apex point 
            faces=[ [0,1,4],[1,2,4],[2,3,4],[3,0,4],              // each triangle side
                    [1,0,3],[2,1,3] ]                         // two triangles for square base
    );

    for (i=[0:1])
        translate([0,0,80*i])
        for (i=[0:3])
            rotate([0,0,i*90])
            translate([0,-40,40])
                window(); 
}

module house4()
{
    tower(1);
    
    translate([0,-25,45]) window();


}

module factory1(w,h)
{
    x = -(w/2-0.5)*100;
    xp = -x;
    y = -(h/2-0.5)*100;
    yp = -y;

    translate([x,y,0])
        wall_edge(1.5);
    
    translate([xp,y,0])
        rotate([0,0,90])
            wall_edge(1.5);

    translate([xp,yp,0])
        rotate([0,0,180])
            wall_edge(1.5);

    translate([x,yp,0])
        rotate([0,0,270])
            wall_edge(1.5);

    translate([0,y,0])
    for (i=side_range(w))
        translate([i*100,0,0])
            wall(1.5);

    translate([0,yp,0])
    for (i=side_range(w))
        translate([i*100,0,0])
            rotate([0,0,180])
                wall(1.5);

    translate([x,0,0])
    for (i=side_range(h))
        translate([0,i*100,0])
            rotate(270)
                wall(1.5);

    translate([xp,0,0])
    for (i=side_range(h))
        translate([0,i*100,0])
            rotate(90)
                wall(1.5);
    for(i=[x,xp], j=[y,yp])
        translate([i,j,0])
            tower(1.75);
}

module factory2(w,h)
{
    x = -(w/2-0.5)*100;
    xp = -x;
    y = -(h/2-0.5)*100;
    yp = -y;

    s = ceil(ln(min(w,h)+3)*2+1);

    difference()
    {
        union()
        {
            for (i=[0:s-1])
                translate([0,0,50*(i+0.5)]) cube([(w-1-i/2)*100,(h-1-i/2)*100,50],center=true);   
        }

        translate([0,0,50*(s-1+0.5)+25]) cube([(w-1-(s-1)/2)*100-50,
                                            (h-1-(s-1)/2)*100-50,
                                            25+$fs], center=true);   
        
    }

    translate([0,0,50*(s-1+0.5)])
    if(min(w,h)>2)
        house2();

}

module factory3(w,h)
{
    range = max( min(min(w,h)/2,4)-1, 0);
    for (s = [0:range])
    {
        wt = w-s;
        ht = h-s;
        xt = -(wt/2-0.5)*100;
        xpt = -xt;
        yt = -(ht/2-0.5)*100;
        ypt = -yt;

        if (wt > 1 && ht > 1)
        {
            translate([xt,yt,0])
                wall_edge(1*s+1.5);
            
            translate([xpt,yt,0])
                rotate([0,0,90])
                    wall_edge(1*s+1.5);

            translate([xpt,ypt,0])
                rotate([0,0,180])
                    wall_edge(1*s+1.5);

            translate([xt,ypt,0])
                rotate([0,0,270])
                    wall_edge(1*s+1.5);

            translate([0,yt,0])
            for (i=side_range(wt))
                translate([i*100,0,0])
                    wall(1*s+1.5);

            translate([0,ypt,0])
            for (i=side_range(wt))
                translate([i*100,0,0])
                    rotate([0,0,180])
                        wall(1*s+1.5);

            translate([xt,0,0])
            for (i=side_range(ht))
                translate([0,i*100,0])
                    rotate(270)
                        wall(1*s+1.5);

            translate([xpt,0,0])
            for (i=side_range(ht))
                translate([0,i*100,0])
                    rotate(90)
                        wall(1*s+1.5);
        }
    }
}

module bighouse1(d)
{
    r = d/2;
    H=ln(d)*150;
    difference()
    {
        cylinder(r=r*100-50, h=H);
        translate([0,0,-$fs]) cylinder(r=r*100-90, h=H+2*$fs);

        for (i=[0:d-1])
        {
            rotate([0,0,180/d*i])
            {
                translate([0,0,25])
                    cube([d*200+$fs,50,50+$fs],center=true);

                translate([0,0,H])
                    cube([d*200+$fs,50,50],center=true);
            }

            rotate([0,0,180/d*(i+0.5)])
            {
                translate([0,0,120])
                    cube([d*200+$fs,50,50],center=true);
            }


        }
    }
}

// factory1(8,8);
// factory2(8,4);
// factory3(16,8);
// bighouse1(4);
// house1();
// house2();
// house3();
// house4();
//bighouse1(8);
//window();
//house4();
//square([800,400],center=true);
