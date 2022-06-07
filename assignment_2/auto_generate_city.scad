include <building_block.scad>

// should not change
block_d = 100;
// should not change
base_h = 50;

isSeperate = true;

isRandomSeed = true;
// if the seed is not random you can animate

city_w = 25;
city_h = 25;
gen_seed = round(rands(0,100000,1)[0]*2+1);

input_seed = 136808;
// city_seed = 140289;
echo(str("The seed is = ", city_seed()));

function city_seed() = isRandomSeed ? gen_seed : input_seed;
// function get_new_seed( s ) = rands(0, 98041, 1)[0] * 2 + 1;
function get_new_seed( s ) = rands(0, 98041, 1, s)[0];
// function get_new_seed( s ) =  ( s * 98041 * 2 + 1 ) % 1403377;
function power_two_less_than( value ) = pow(2, floor( ln(value) / ln(2) ));

// unit to real
function conv( unit ) = unit * block_d;
function convpos( unit ) = isSeperate ? unit * block_d * (1 + 0.5/($t*10+1)) : conv(unit);

function rands10( s ) = rands(0,10,1,seed=s)[0];


module city( w, h, seed )
{
    seed = get_new_seed(seed);
    city_sec( 0, 0, w, h, seed);
}

module city_sec( x, y, x_p, y_p, seed )
{
    seed1 = get_new_seed(seed);
    seed2 = get_new_seed(seed1);
    seed3 = get_new_seed(seed2);

    w = x_p-x;
    h = y_p-y;
    if (w > 0 && h > 0)
    {
        max_sq_d = min( power_two_less_than(w), power_two_less_than(h) );

        // square module
        sq_sec( x, y, max_sq_d, seed1 );

        // lower module
        city_sec( x, y + max_sq_d, x + max_sq_d, y_p, seed2 );

        // righter module
        city_sec( x + max_sq_d, y, x_p, y_p, seed3 );

    }

}

module sq_sec( x, y, d, seed )
{
    rand = rands10(seed);

    seed1 = get_new_seed(seed);
    seed2 = get_new_seed(seed1);

    rec_p = 15/pow(1.3,d);
    bar_p = (10-rec_p)/2;

    if ( rand > 10-rec_p-bar_p && d>1 )
    {
        h_sec( x, y, d, seed1 );
        h_sec( x, y+d/2, d, seed2 );
    }
    else if ( rand > rec_p && d>1 )
    {
        v_sec( x, y, d, seed1 );
        v_sec( x+d/2, y, d, seed2 );
    }
    else
    {
        seed1 = get_new_seed(seed);
        sq_block( x, y, d, seed1);
    }

}

module h_sec(x,y,d,seed)
{
    rand = rands10(seed);
    seed1 = get_new_seed(seed);
    seed2 = get_new_seed(seed1);
    seed3 = get_new_seed(seed2);

    rec_p = 5.0/pow(1.1,d);

    if (rand > 10-rec_p && d>2 )
    {
        h_block(x,y,d,seed3);
    }
    else
    {
        sq_sec(x,y,d/2,seed1);
        sq_sec(x+d/2,y,d/2,seed2);
    }
}

module v_sec(x,y,d,seed)
{
    rand = rands10(seed);
    seed1 = get_new_seed(seed);
    seed2 = get_new_seed(seed1);
    seed3 = get_new_seed(seed2);

    rec_p = 5.0/pow(1.1,d);

    if (rand > 10-rec_p && d>2 )
    {
        v_block(x,y,d,seed3);
    }
    else
    {
        sq_sec(x,y,d/2,seed1);
        sq_sec(x,y+d/2,d/2,seed2);
    }
}

module sq_block(x,y,d,seed)
{
    seed1 = get_new_seed(seed);
    randh = rands10(seed1);

    base(x,y,x+d,y+d);
    translate([convpos(x) + block_d*d/2, convpos(y)+block_d*d/2,base_h])
    {
        if (d==1)
        {

            if (randh > 8.5)
                {house1();}
            else if (randh > 7)
                {house2();}
            else if (randh > 5.5)
                {house3();}
            else if (randh > 4)
                {house4();}
            

        }
        else
        {

            if (randh > 7.5)
            {
                bighouse1(d);
            }
            else if (randh > 5)
            {
                factory1(d,d);
            }
            else if (randh > 2.5)
            {
                factory2(d,d);
            }
            else{
                factory3(d,d);
            }
        }
    }
}

module h_block(x,y,d,seed)
{
    seed1 = get_new_seed(seed);
    rand = rands10(seed1);
    w = d;
    h = d/2;

    translate([convpos(x) + block_d*w/2, convpos(y)+block_d*h/2,base_h])
    {

        if (rand >0.666)
            factory1(w,h);
        else if (rand > 0.333)
            factory2(w,h);
        else
            factory3(w,h);

    }
        
    base( x, y, x + w, y + h );
    
}

module v_block(x,y,d,seed)
{
    seed1 = get_new_seed(seed);
    rand = rands10(seed1);
    w = d/2;
    h = d;

    translate([convpos(x) + block_d*w/2, convpos(y)+block_d*h/2,base_h])
    {

        if (rand >0.666)
            factory1(w,h);
        else if (rand > 0.333)
            factory2(w,h);
        else
            factory3(w,h);

    }
        
    base( x, y, x + w, y + h );
}

module base( x, y, x_p, y_p )
{
    w = x_p - x;
    h = y_p - y;

    joint_r = base_h * 0.4;
    joint_l = block_d / 4;
    translate([convpos(x), convpos(y),0]){
        difference()
        {
            cube( [ conv(w), conv(h), base_h ] );

            for ( i=[0:w-1], j=[0:h-1] )
            {
                if (y != 0)
                    translate([ conv(i+0.5), 0-$fs, base_h/2 ]) rotate([-90,0,0]) cylinder(r=joint_r, h=joint_l);
                if (x != 0)
                    translate([ 0-$fs, conv(j+0.5), base_h/2 ]) rotate([0,90,0]) cylinder(r=joint_r, h=joint_l);
            }
        }


        for ( i=[0:w-1], j=[0:h-1] )
        {
            if (y_p != city_h)
                translate([ conv(i+0.5), conv(h), base_h/2 ]) rotate([-90,0,0]) cylinder(r=joint_r*0.85, h=joint_l);
            if (x_p != city_w)
                translate([ conv(w), conv(j+0.5), base_h/2 ]) rotate([0,90,0]) cylinder(r=joint_r*0.85, h=joint_l);
        }
    }
}
color("grey") city( city_w, city_h, city_seed() );
// cylinder(r=10,h=20);
// echo(str("Variable = ", rands(0,10,1,seed=111)));
// base(0,0,4,4);
// sq();