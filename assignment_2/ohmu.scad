elevate = 1/3;

function ngon(num, r) = 
    [for (i=[0:4], a=i*360/num) [ r*cos(a), r*sin(a) ]];

function partial_hexagon(r) =
    [for (i=[0:6], a=i*360/8) [ r*cos(a), r*sin(a) ]];

module base ( r )
{
    rotate([0,0,-45])
    polygon(partial_hexagon(r));

}

module extrude_partial_hexagon( r, R, l )
{
    scale_factor = R/r;
    linear_extrude(height = l, scale = scale_factor )
    base( r );

}

module _segment_joint_female( d, scale_factor, h )
{
    linear_extrude(height=h, scale=[1/scale_factor,1/scale_factor])
        circle(d=d*scale_factor);
    translate([0,0, h+d/2])
        cube([d * 1.5 , d*1.5, d], center=true);
}

module _segment_joint_male( d, l )
{
    cylinder(r=d/4, h=l+d/2);

    translate([0,0,l+d/2])
        sphere(d=d);
}

module segment( r, R, l )
{
    shell_p = 1/8;
    inside_p = 1-shell_p;
    shell_deep_p = 1/3;

    shell_thick = (R-r) * shell_deep_p;
    shell_deep = l*shell_deep_p;
    
    seg_gap = shell_deep*1/2;

    inner_section_r = ( r + shell_thick ) * inside_p;

    // joint female
    hole_d = 50;


    difference()
    {   
        // core
        translate( [0,0,-l*shell_deep_p] )
            difference()
            {
                extrude_partial_hexagon(r, R, l);

                translate([0,0,-$fs]) 
                    extrude_partial_hexagon(r*inside_p, inner_section_r, shell_deep +$fs);

                translate([-R,-R, l*shell_deep_p - R])
                    cube([2*R,R,R]);
            
            }
        
        
        _segment_joint_female(d=hole_d, scale_factor=3, h=shell_thick);

    }

    male_joint_offset = (R-r)*seg_gap/l;
    translate([0,male_joint_offset, l - shell_deep])
        _segment_joint_male(hole_d, seg_gap);
    
        

}

module seg_rec(r,R,l, divY)
{
    translate([0,0,divY])
        segment(R=R, r=r, l=l);
    
    if (r > 100)
    {
        seg_rec(R=r, r=r-100, l=l, divY=divY-l);
    }
}

module model()
{
    seg_rec(r=500,R=600, l = 1000, divY=0);
}

model();
// segment( r=100, R=200, l=300);
// _segment_joint_male(d=25, l =100);
//_segment_joint_female(25,5,20);

