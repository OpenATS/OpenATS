/*
欢迎使用OpenATS天线追踪系统（Open Antenna Tracking System）



需要的硬件：
  电脑一台，用来运行WXtrack卫星跟踪软件。
  Arduino，山寨或者正版都可以，正版可以购买Genuino，分Nano和Mega等版本，建议Mega，拥有更好的处理能力和内存，并且以后扩展好。国内正版版本大概140元左右。Nano大概80左右，山寨30-50左右。山寨也正常，我的实验平台就是在山寨上运行的
  步进电机2个，最好用57步进电机，扭距选择大一些的，2.0Nm以上最好，当然还可以通过减速机来增大扭矩。
  步进电机驱动器2个，分别控制步进电机的。1套驱动器配1个步进电机，57步进套装在200-300左右，当然可以淘二手的。
  24V开关电源一个，工业上用的即可，淘宝很多，具体功率要根据你两个步进电机的功率之和来决定。如果选57大功率的步进电机的话，就要选择最少10A以上的电源
。
  再一个不起眼但很重要的，XY轴可以旋转的支架。这个不好购买，可以自己制作，成本不太好估计。
  如果不算支架的话，仅需要几百块钱就可以完成所实现的功能，加上支架成本，不到2k应该足够了。具体看你用没用减速机还有支架配件等等，有的配件，我们身边都有。


搭建过程：

所有程序和代码，已经分享到我的百度网盘，地址：http://pan.baidu.com/s/1i5jR1I1


1，把Arduino驱动程序arduino-1.6.8-windows.exe和虚拟串口软件驱动CH341SER.EXE下载安装，把AccelStepper库下载解压缩，解压缩后的整个AccelStepper文件夹拷贝到安装后Arduino程序的libraries文件夹下
   将追踪天线代码上传到Arduino板子中。

2，搭建硬件，具体硬件接线图，可以看图片。简单就是用Arduino的PWM接口，来发送脉冲。AZ方位角的控制接口是3，dir方向控制是5，EL仰角的控制接口是6，dir方向控制是9。采用共阴接法或者共阳接法都可以。
   
   AZ方位角的CLK+(有的驱动器上也叫PUL+)接到Arduino的数字接口3，CW+(有的也叫DIR+）接到数字口5，EL仰角的CLK+接到6，CW+接到9上
   AZ方位角驱动器上的CLK-，CW-接到一起，再接到arduino的一个GND接口。EL仰角驱动器的CLK-，CW-接到一起，再接到Arduino的另一个GND接口。
   剩余的EN+（有的也叫ENA+）还有EN-我们都不接。

3，步进电机与驱动器的接线，步进电机比如常见的2相四线电机，分A+,A-,B+,B-接线，在驱动器上都有标识，自行按照自己步进电机的线颜色定义接上即可。驱动器上的电源接开关电源的输出即可，开关电源给两个驱动器供电24V电压。
   注意正负极。



使用说明：

（一）自动追踪

  自动追踪，电脑下载运行WXtrack软件，更新星历，星历更新会提示联不上网或者更新失败，我们手动去下载星历文件放到程序的文件夹下面即可。
  电脑运行WXtrack软件，找到顶部菜单，Tracker--Options..然后打开设置天线追踪器界面，然后我们点击选中EasyComm I COM port，
  下面的窗口，pass start可以设置1，意思是1分钟之前就唤醒天线。Pass end就是追踪完成后发送什么指令，我们输入我们天线的复位指令：“S”，
  Parking意思是停止天线，我们选中Park antenna，其余不用改。再切换到Port setup菜单，选中你的arduino对应的COM口，比如COM4，速率选中程序中设置的9600
  然后是进入Tracker-specific options菜单，将EasyComm precision（精确度控制）设置为2，这样程序输出的串口数据会保留小数点2位，否则是整数，整个系统精确度很低。
  电脑软件我们就设置完成了。调试好角度后，就可以自动追踪啦！

（二）手动控制

如果手动控制，请打开Arduino的IDE自带的串口监视器或者用别的可以发送串口数据的软件（各种串口调试工具都可以），用来做手动控制

设置你的Arduino对应的串口编号，然后设置波特率9600。
在发送角度命令时，按照如下格式：方位角 仰角（AZ EL），中间为空格，数据为浮点数或者整数，比如发送：20 40，则天线方位角转到20度，仰角到40度。
等同于：20.0 40.0  可以输入小数点儿，细分大的话，程序的识别能力可以精确到0.1度。
如果仅仅输入一位数字，则只调动方位角，仰角便为0.比如：20则天线默认将方位角转到20度，仰角为0度。输入R或者0都可以将天线复位，即方位角仰角都为0度。

以下功能没测试，不能正常运行请见谅，待以后测试更新。
校准功能
在长时间运转后，天线可能会由于电机的丢步，造成的精度不准，这时候可以输入命令来进行校准。为了简单明了来调整天线，AZ跟EL分开调整。
调整方位角的命令为：AA23 就是大写字母A后面紧跟一个数字，可以是浮点数也可以是整数。就是把天线当前的角度调整为23度。
调整仰角的命令便为：EE56 原理同上。




后续：

本人水平有限，非程序员，只是个普通工人，工作繁忙。Arduino的代码也没有仔细学，只是为了达到我的目的而去补的知识
代码如果不好，请见谅。公开代码就是让大家都能自己制作自己的追踪系统，不必花那么多钱买商业的成品。代码中已经注释的非常详细

希望更多的人完善我们的天线，来将我们的天线越做越好。在此祝福大家在玩无线电的过程中享受乐趣。

73

*/
/////////////////////////////////////////////////////////////////////////////////////////////

#include <AccelStepper.h>

#define ACCELERATION_1 400.0                      //AZ步进电机加速度参数
#define ACCELERATION_2 400.0                      //EL步进电机加速度参数
#define MAXSPEED_1 700.0                          //AZ步进电机最大速度
#define MAXSPEED_2 700.0                          //EL步进电机最大速度
#define STEPPER1_STEP_PIN 3                       //stepper1 AZ步进电机脉冲针脚
#define STEPPER1_DIR_PIN  5                       //AZ方向控制针脚
#define STEPPER2_STEP_PIN 6                       //stepper2 EL步进电机脉冲针脚
#define STEPPER2_DIR_PIN  9                       //EL方向控制针脚


#define AZFACTOR 44.444                           //AZ脉冲角度比         这是临时假设的值，没有细调
#define ELFACTOR 44.444                           //EL脉冲角度比         这是临时假设的值，没有细调

// 上面两个角度比函数的数值比较重要，不是固定的，而是根据你的天线的脉冲数和角度值的比值得来的
// 假如你的驱动器设置步进电机的细分为1的话，两相步进电机则为200步一圈，而一圈AZ为360度，则得出比例为200/360=0.55555...意思为0.5555个脉冲则旋转1度。
// 如果你加了1：5的减速机，则为200*5=1000步/圈，比例系数则为1000/360=2.77777，如果设置了8细分，则为：200*8*5/360=22.222222
// 这只是理论值，但是安装好天线后，还要进行调试，因为机械误差等等，不一定会正好旋转一圈，具体数值还要根据你所用的天线机械装置来确定。
// 尽量将细分设置大些，比如8细分，16细分。如果普通驱动器的话，不建议设置更大的细分，那样会造成步进电机的丢步。
// 丢步为步进电机的比较严重的问题，因为我们不能读出步进电机的角度值，所以步进电机丢步后，我们却不知道，会造成以后所有的跟踪都不准
// 速度设置慢些，加速度值不要太大，天线重量要平衡，这样可以防止步进电机的丢步。还要定时校准一下天线，如果有条件，就上伺服电机，当然上了伺服，这些代码则失效了。

float X;
float Y;
float W;
float M;
float gotoangle_x;                                 //步进电机要走的方位角脉冲
float gotoangle_y;                                 //步进电机要走的仰角脉冲
float angle_x_tmp = 0;                             //方位角缓存参数用来确定AZ角差值

String usbdata = "";                               //定义一个空字符串
boolean _angle_1 = false;                          //判断方位角是否符合要求的布尔值
boolean _angle_2 = false;                          //判断仰角是否符合要求的布尔值

int angle_x  = 0;                                  //差值比较函数，当卫星横跨了0度，则激活某个条件，再后续的追踪，都走另一种算法。条件有（-1、0、1）两种情况，0为默认

/////////////////////////-以下为声明步进电机针脚-/////////////////////////////////////////
AccelStepper stepper1(AccelStepper::DRIVER, STEPPER1_STEP_PIN, STEPPER1_DIR_PIN);
AccelStepper stepper2(AccelStepper::DRIVER, STEPPER2_STEP_PIN, STEPPER2_DIR_PIN);
//////////////////////////////////////////////////////////////////////////////////////////

void setup()
{  
   



   float angle_x_tmp = 0;                         //方位角缓存参数用来确定AZ角差值，用来缓存上一次的方位角度，当前角度和上次角度做对比来识别卫星是否是横跨了0度方位角
                                                  //如果卫星从北方（0度）线跨过时，差值会明显很大，分两种情况，然后根据不同情况来使天线倒转，防止天线自动再转一圈来追踪，造成信号有盲区
                                                  //比如一个卫星，从0度飞到了359度时，检测到差值符合条件，则天线直接倒转到-1度，防止天线再正传1圈到359度，要不然这段时间就是盲区了
                                                  //反过来如果从359度到0度时，就会自动转到360度，然后接着361度、362度继续转防止天线再回转一圈。追踪完成复位后才不影响后续追踪，可以手动复位发送“S”



   //stepper1电机为方位角AZ
   //stepper2电机为仰角  EL



   stepper1.setMinPulseWidth(10);                 //脉冲宽度（可以控制速度）
   stepper2.setMinPulseWidth(10);                 //脉冲宽度（可以控制速度）

   stepper1.setAcceleration(ACCELERATION_1);
   stepper1.setMaxSpeed(MAXSPEED_1);

   stepper2.setAcceleration(ACCELERATION_2);
   stepper2.setMaxSpeed(MAXSPEED_2);

 //stepper1.setSpeed(1800.0);
 //stepper2.setSpeed(1800.0);
  
   Serial.begin(9600);
   Serial.println("The Open Auto Tracking System is online !");
   Serial.println("Use 'XX YY' to control the system");
   Serial.println("XX is AZ, YY is EL");
   Serial.println("Like this : 35.6 45.0");
   Serial.println("Press 'S' or '0' to home"); 
   Serial.println("Press 'AA52.4' will calibration AZ to 52.4");
   Serial.println("Press 'EE30.6' will calibration EL to 30.6");
   Serial.println("Good Luck !!!");
}
//////////////////////////////////////////////////////////////////////////////////////////

void loop()                                        //主程序循环
{
//下面为初始化各种变量，请勿移动位置

String usbdata = "";
String tmp_a = "";
String tmp_b = "";
String tmp_c = "";
String tmp_d = "";
String tmp_e = "";
String tmp_f = "";
char a[10]={0};
char b[10]={0};
char c[10]={0};
char d[10]={0};
char e[10]={0};
char f[10]={0};

//下面为接收串口数据，将串口字符串数据，传递给usbdata变量

   while (Serial.available() > 0) 
   {
    usbdata += char(Serial.read());
    delay(2);
   }                                                                    

//自动追踪下usb串口的数据格式例子：

//Easycomm I COM port模式的数据为：（推荐注册后的模式）
//AZ23.0 EL56.0 UP145978126 DN432115748\#13
//AZ84.46 EL19.41\#13（软件注册后）

//KVH模式下的串口数据为：
//AZ,239\#13\#10       方位角23.9度                   
//EL,241\#13\#10       仰角  24.1度

//Easycomm模式方位角跟仰角在同一个串口命令中，这样两电机可以同时旋转，而KVH模式则是分开发送2个角度数据
//虽然都可以完成任务，但是自动追踪的过程中，KVH模式两个电机会不同时旋转，这样对精度就会造成较大的影响
//并且追踪系统运行的过程中不太美观，但是一个比较重要的一点儿是，由于WXtrack需要注册才能更好的输出EasyComm
//更精准的浮点数(即将推出的测试版将免费),而网上有WXtrack 3.8.28版本的破解注册机，需要破解后也可以完成，但破解对于原作者的贡献实在是愧疚
//而这个天线系统不能做那样的宣传，请低调使用，如果你条件允许，请支持原作者的付出。如果你有需要请自由切换
//仅需要把两个模式的代码注释掉一种即可。




////////////////////////////////////////////////////////////////////////////////
//-------------------下面代码为追踪器选择的是Easycomm模式---------------------//
////////////////////////////////////////////////////////////////////////////////
//软件注册后的串口数据格式为：AZ84.46 EL19.41\#13

 if (usbdata.length() > 15)
   {

      const char *usb=usbdata.c_str();                           //将string数据转成字符串数据
      //Serial.println(usb);
      sscanf(usb, "AZ%s EL%s\#13", &a, &b);                      //进行匹配字符串数组
      tmp_a = a;                                                
      tmp_b = b;                                                 
                                          
      X = tmp_a.toFloat();                                       //进行浮点数转换
      Y = tmp_b.toFloat();                                      
      //Serial.println(X); 
      //Serial.println(Y); 
      if( X-angle_x_tmp < -350)                                  //此判断式作用是检测是否是从359度过渡到0度，不至于方位角倒转一圈
          {
             angle_x = 1;      
          }
      else if( X-angle_x_tmp > 350 )                             //此判断式作用是检测是否是从0度过渡到359度，不至于方位角倒转一圈
          {
             angle_x = -1;
          }

      angle_x_tmp = X;                                           //重新赋值新的比较角度
      gotoangle_x = X;
      gotoangle_y = Y;
      _angle_1 = true;   
      _angle_2 = true;

   }
 else if(usbdata.startsWith("S"))                                //如果数据开头是停止命令“S”
   {
      gotoangle_x = 0;
      gotoangle_y = 0;
      angle_x = 0;
      _angle_1 = true;
      _angle_2 = true;
   } 
   
 else
   {                                                             //否则进行手动控制
  
    if(usbdata.startsWith("AA"))                                 //检测是否为校准方位角命令
      { 
      const char *usb=usbdata.c_str();                          
      sscanf(usb, "AA%s",&e);
      tmp_e = e;
      W = tmp_e.toFloat();
      stepper1.setCurrentPosition(W);
    
      }

    else if(usbdata.startsWith("EE"))                            //检测是否为校准仰角命令
      { 
      const char *usb=usbdata.c_str();                          
      sscanf(usb, "EE%s",&f);
      tmp_f = f;
      M = tmp_f.toFloat();
      stepper2.setCurrentPosition(M);
    
      }

    else if( usbdata.length() > 0)                                
     {
      const char *usb=usbdata.c_str();                           //将string数据提取仰角数据
      sscanf(usb, "%s %s",&c, &d);
      tmp_c = c;
      tmp_d = d;
      X = tmp_c.toFloat();
      Y = tmp_d.toFloat();
      if (X < 400 && X > -400)                                   //手动控制时，进行角度限制。
      gotoangle_x = X;
      if (Y > -20 && Y < 95)       
      gotoangle_y = Y;
      angle_x = 0;
      _angle_1 = true;
      _angle_2 = true;
    }

   }
////////////////////////////////////////////////////////////////////////////////
//------------------Easycomm追踪模式的代码结束--------------------------------//
////////////////////////////////////////////////////////////////////////////////


/*

/////////////////////////////////////////////////////////////////////////////////
//-------------------下面代码为追踪器选择的是KVH模式---------------------------//
/////////////////////////////////////////////////////////////////////////////////
//KVH模式下的串口数据为：
//AZ,239\#13\#10       方位角23.9度                   
//EL,241\#13\#10       仰角 24.1度


//下面为判断usbdata第一个字符，来分辨是自动控制还是手动控制

if(usbdata.startsWith("A")||usbdata.startsWith("E"))
{
    if (usbdata.startsWith("A"))
       {
             const char *usb=usbdata.c_str();                   //将string数据提取方位角数据
             sscanf(usb, "AZ,%s\#13\#10", &a);                  //进行匹配字符串数组
             tmp_a = a;
             X = tmp_a.toFloat()/10;                            //浮点数转换成角度值  
                                                                //假如你的机械结构是反转，则加上“-”号，即X=-(tmp_a.toFloat()/10)即可，相应的Y也是同样的道理

        if( X-angle_x_tmp < -350)                               //此判断式作用是检测是否是从359度过渡到0度，不至于方位角倒转一圈
          {
             angle_x = 1;      
          }
        else if( X-angle_x_tmp > 350 )                          //此判断式作用是检测是否是从0度过渡到359度，不至于方位角倒转一圈
          {
             angle_x = -1;
          }

      angle_x_tmp = X;                                          //重新赋值新的比较角度
      gotoangle_x = X;
      _angle_1 = true;   
      
      //Serial.println(X);

      }
   else  // if (usbdata.startsWith("E"))
      {
      const char *usb=usbdata.c_str();                          //将string数据提取仰角数据
      sscanf(usb, "EL,%s\#13\#10", &b);
      tmp_b = b;     
      Y = tmp_b.toFloat()/10;      
      gotoangle_y = Y; 
      _angle_2 = true;
      //Serial.println(Y);
      }
}
 else if(usbdata.startsWith("S"))                               //如果数据开头是停止命令“S”
      {
      gotoangle_x = 0;
      gotoangle_y = 0;
      angle_x = 0;
      _angle_1 = true;
      _angle_2 = true;
      } 
   
 else
   {                                                            //否则进行手动控制
  
    if(usbdata.startsWith("AA"))                                //检测是否为校准方位角命令
      { 
      const char *usb=usbdata.c_str();                          
      sscanf(usb, "AA%s",&e);
      tmp_e = e;
      W = tmp_e.toFloat();
      stepper1.setCurrentPosition(W);
    
      }

    else if(usbdata.startsWith("EE"))                           //检测是否为校准仰角命令
      { 
      const char *usb=usbdata.c_str();                          
      sscanf(usb, "EE%s",&f);
      tmp_f = f;
      M = tmp_f.toFloat();
      stepper2.setCurrentPosition(M);
    
      }

    else if( usbdata.length() > 0)                                
      {
      const char *usb=usbdata.c_str();                          //将string数据提取仰角数据
      sscanf(usb, "%s %s",&c, &d);
      tmp_c = c;
      tmp_d = d;
      X = tmp_c.toFloat();
      Y = tmp_d.toFloat();
      if (X < 400 && X > -400)                                  //手动控制时，进行角度限制。
      gotoangle_x = X;
      if (Y > -20 && Y < 95)       
      gotoangle_y = Y;
      angle_x = 0;
      _angle_1 = true;
      _angle_2 = true;
      }

       
}
/////////////////////////////////////////////////////////////////////////////////
//----------------------------KVH追踪模式代码结束------------------------------//
/////////////////////////////////////////////////////////////////////////////////

*/



//方位角旋转角度

   if (_angle_1)     
     {
        if( angle_x > 0 )                                       //判断是否是从0度过渡到359.x度
        {       
        stepper1.moveTo((360+gotoangle_x) * AZFACTOR);   

        }
        else if( angle_x < 0 )                                  //方位角从359.x过渡到0度   
        {    
        stepper1.moveTo((gotoangle_x-360) * AZFACTOR);  

        }
        else                                                  
        {       
        stepper1.moveTo(gotoangle_x * AZFACTOR);   
        }
                                                 

     }
//仰角旋转角度
   
   if (_angle_2)   
       stepper2.moveTo(gotoangle_y * ELFACTOR);   
       stepper1.run();
       stepper2.run();
    
}
//整个程序结束

