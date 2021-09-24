/*
欢迎使用OpenATS自动追踪天线（Open Auto Tracking System）

(WXtrack控制软件）

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

所有程序和代码，已经分享到我的百度网盘，链接:https://pan.baidu.com/s/1sliIPJJ 密码:pdx2 链接:https://pan.baidu.com/s/1dZcN3S 密码:7sxy https://pan.baidu.com/s/1nvtFfTJ 解压密码：1024


1，把Arduino驱动程序arduino-1.6.8-windows.exe和虚拟串口驱动CH341SER.EXE下载安装，把AccelStepper库下载解压缩，解压缩后的整个AccelStepper文件夹拷贝到安装后Arduino程序的libraries文件夹下
   将追踪天线代码上传到Arduino板子中。

2，搭建硬件，具体硬件接线图，可以看图片。简单就是用Arduino的PWM接口，来发送脉冲。AZ方位角的控制接口是5，dir方向控制是6，EL仰角的控制接口是9，dir方向控制是10。采用共阴接法或者共阳接法都可以。
   
   AZ方位角的CLK+(有的驱动器上也叫PUL+)接到Arduino的数字接口5，CW+(有的也叫DIR+）接到数字口6，EL仰角的CLK+接到9，CW+接到10上
   AZ方位角驱动器上的CLK-，CW-接到一起，再接到arduino的一个GND接口。EL仰角驱动器的CLK-，CW-接到一起，再接到Arduino的另一个GND接口。
   剩余的EN+（有的也叫ENA+）还有EN-我们都不接。（请注意与之前发布的3、5、6、9接口不同，预留出2、3接口目的以后研发更多功能采用外部中断。）

3，步进电机与驱动器的接线，步进电机比如常见的2相四线电机，分A+,A-,B+,B-接线，在驱动器上都有标识，自行按照自己步进电机的线颜色定义接上即可。驱动器上的电源接开关电源的输出即可，开关电源给两个驱动器供电24V电压。
   注意正负极。



使用说明：

（一）自动追踪

  自动追踪，电脑下载运行WXtrack软件，更新星历，星历更新会提示联不上网或者更新失败，我们手动去下载星历文件放到程序的文件夹下面即可。
  电脑运行WXtrack软件，找到顶部菜单，Tracker--Options..然后打开设置天线追踪器界面，然后我们点击选中CX6DD WispDDE Client，
  下面的窗口，pass start,发送命令设置天线的唤醒命令“W”Pass end就是追踪完成后发送什么指令，我们输入我们天线的复位指令：“S”，
  Parking意思是停止天线，我们选中Park antenna，其余不用改。再切换到Port setup菜单，选中你的arduino对应的COM口，比如COM4，速率选中程序中设置的19200
  然后是进入Tracker-specific options菜单，将EasyComm precision（精确度控制）设置为2，这样程序输出的串口数据会保留小数点2位，否则是整数，整个系统精确度很低。如果使用FastWXtrack可以忽略
  设置好自己天线的经纬度坐标，海拔等信息，打开FastWXtrack自动识别，选中卫星，再打开天线的控制客户端，OpenATS DDE Client，所有选项全部选中。即可开始自动控制天线。
  电脑软件我们就设置完成了。调试好角度后，就可以自动追踪啦！

（二）手动控制

如果手动控制，请打开Arduino的IDE自带的串口监视器或者用别的可以发送串口数据的软件（各种串口调试工具都可以），用来做手动控制。设置你的Arduino对应的串口编号，然后设置波特率38400。
默认情况，打开串口通信软件后，天线自动通电。
如果没有的话，会提示天线电源关闭，发送一个W等待两秒就可以了。

在发送角度命令时，按照如下格式：方位角 仰角（AZ EL），中间为空格，数据为浮点数或者整数，比如发送：20 40，则天线方位角转到20度，仰角到40度。
等同于：20.0 40.0  可以输入小数点儿，细分大的话，程序的识别能力可以精确到0.1度。
如果仅仅输入一位数字，则只调动方位角，仰角便为0.比如：20则天线默认将方位角转到20度，仰角为0度。输入R或者0都可以将天线复位，即方位角仰角都为0度。
输入空格、或者任何一个小写字母，天线归零。

*****重要*****
发送S命令，会将天线电源切断，S命令接收后，等待一段时间，等待天线归零后，自动切断电源。此时天线属于睡眠状态。如果想重新使用天线，请发送W命令打开天线的电源。 
如果软件自动追踪，会自动开启电源。
发送L命令，会将天线锁定当前角度，切断电机电源，LNA电源不切断，用于接收固定轨道卫星时节电。发送U命令解除锁定。

校准功能
在长时间运转后，天线可能会由于电机的丢步，造成的精度不准，这时候可以输入命令来进行校准。为了简单明了来调整天线，AZ跟EL分开调整。
调整方位角的命令为：X 23 就是大写字母X后面紧跟一个数字，可以是浮点数也可以是整数。就是把天线的23度调整为0度。
调整仰角的命令便为：Y 56 原理同上。

安装天线时，需要将天线的0度对准地球的正北，而非磁北，所以用指南针安装会出现误差。可以想别的办法，比如借助谷歌地球，获得正北方向的参照物，然后天线直线对准。
或者参照太阳在正南时，所照射的投影为正北方向。也可以使用高精度的GPS来测试出同一条经线。通过WXtrack的Tracker--Test选项，选中sun（太阳），进行追踪太阳，看天线的天馈影子是否在
天线正中心，阳光聚焦点是否在天线天馈的地方。



*/
/////////////////////////////////////////////////////////////////////////////////////////////
//#include<pt.h>                                  //超级多线程库，本程序暂时未使用，以后研发新功能将使用
#include<AccelStepper.h>

#define ACCELERATION_1 400.0                      //AZ步进电机加速度参数
#define ACCELERATION_2 400.0                      //EL步进电机加速度参数
#define MAXSPEED_1 700.0                         //AZ步进电机最大速度
#define MAXSPEED_2 700.0                         //EL步进电机最大速度
#define STEPPER1_STEP_PIN 5                       //stepper1 AZ步进电机脉冲针脚
#define STEPPER1_DIR_PIN  6                       //AZ方向控制针脚
#define STEPPER2_STEP_PIN 9                       //stepper2 EL步进电机脉冲针脚
#define STEPPER2_DIR_PIN 10                       //EL方向控制针脚


#define AZFACTOR 26.666                           //AZ脉冲角度比         这是临时假设的值，没有细调
#define ELFACTOR 26.666                           //EL脉冲角度比         这是临时假设的值，没有细调

// 上面两个角度比函数的数值比较重要，不是固定的，而是根据你的天线的脉冲数和角度值的比值得来的
// 假如你的驱动器设置步进电机的细分为1的话，两相步进电机则为200步一圈，而一圈AZ为360度，则得出比例为200/360=0.55555...意思为0.5555个脉冲则旋转1度。
// 如果你加了1：5的减速机，则为200*5=1000步/圈，比例系数则为1000/360=2.77777，如果设置了8细分，则为：200*8*5/360=22.222222
// 这只是理论值，但是安装好天线后，还要进行调试，因为机械误差等等，不一定会正好旋转一圈，具体数值还要根据你所用的天线机械装置来确定。
// 尽量将细分设置大些，比如8细分，16细分。如果普通驱动器的话，不建议设置更大的细分，那样会造成追踪速度降低。
// 丢步为步进电机的比较严重的问题，因为我们不能读出步进电机的角度值，所以步进电机丢步后，我们却不知道，会造成以后所有的跟踪都不准
// 速度设置慢些，加速度值不要太大，天线重量要平衡，这样可以防止步进电机的丢步。还要定时校准一下天线，如果有条件，就上伺服电机，当然上了伺服，这些代码则失效了。
//本代码中的 26.666计算由来：2相步进电机，200*6倍减速比*8细分/360=26.666。有不懂的，或者不明白的，请加微信群聊。


int power;
int looptime = 1;
int angle_x   = 0;                                 //差值比较函数，当卫星横跨了0度，则激活某个条件，再后续的追踪，都走另一种算法。条件有（-1、0、1）两种情况，0为默认
int power1Pin = 4;
int power2Pin = 8;
int power_tmp = 0;                                 //0为正常状态，-1为锁定状态，1为断电状态

float X;
float Y;
float M;
float N;
float gotoangle_x;                                 //步进电机要走的方位角脉冲
float gotoangle_y;                                 //步进电机要走的仰角脉冲
float angle_x_tmp = 0;                             //方位角缓存参数用来确定AZ角差值

String usbdata = "";                               //定义一个空字符串
boolean _angle_1 = false;                          //判断方位角是否符合要求的布尔值
boolean _angle_2 = false;                          //判断仰角是否符合要求的布尔值

//unsigned long stoptime;
//unsigned long starttime;                         //计算代码循环一次时间函数
unsigned long accumulate;                          //循环时间累加函数
//static struct pt pt1,pt2;
/////////////////////////-以下为声明步进电机针脚-/////////////////////////////////////////
AccelStepper stepper1(AccelStepper::DRIVER, STEPPER1_STEP_PIN, STEPPER1_DIR_PIN);
AccelStepper stepper2(AccelStepper::DRIVER, STEPPER2_STEP_PIN, STEPPER2_DIR_PIN);
//////////////////////////////////////////////////////////////////////////////////////////

void setup()
{  
   
   pinMode(power1Pin, OUTPUT);
   pinMode(power2Pin, OUTPUT);
   //PT_INIT(&pt1);                                  
   //attachInterrupt(2, power1, LOW);             //外部中断定义
   //attachInterrupt(3, power2, LOW);
   float angle_x_tmp = 0;                         //方位角缓存参数用来确定AZ角差值，用来缓存上一次的方位角度，当前角度和上次角度做对比来识别卫星是否是横跨了0度方位角
                                                  //如果卫星从北方（0度）线跨过时，差值会明显很大，分两种情况，然后根据不同情况来使天线倒转，防止天线自动再转一圈来追踪，造成信号有盲区
                                                  //比如一个卫星，从0度飞到了359度时，检测到差值符合条件，则天线直接倒转到-1度，防止天线再正传1圈到359度，要不然这段时间就是盲区了
                                                  //反过来如果从359度到0度时，就会自动转到360度，然后接着361度、362度继续转防止天线再回转一圈。追踪完成复位后才不影响后续追踪，可以手动复位发送“S”

  
   //stepper1电机为方位角AZ
   //stepper2电机为仰角  EL


   
   stepper1.setMinPulseWidth(10);                  //脉冲宽度（可以控制速度）
   stepper2.setMinPulseWidth(10);                  //脉冲宽度（可以控制速度）

   
   stepper1.setMaxSpeed(MAXSPEED_1);
   stepper2.setMaxSpeed(MAXSPEED_2);
   
   stepper1.setAcceleration(ACCELERATION_1);
   stepper2.setAcceleration(ACCELERATION_2);
   

   //stepper1.setSpeed(8000.0);
   //stepper2.setSpeed(8000.0);
  
   Serial.begin(9600);
   
   Serial.println("++++++++++++++++++++++++++++++++++++++++++++++++++++");
   Serial.println("+    The Open Auto Tracking System is online !     +");
   Serial.println("+                   (WXtrack)                      +");
   Serial.println("+        Use 'XX YY' to control the antenna        +");
   Serial.println("+               XX is AZ, YY is EL                 +");
   Serial.println("+           Like this : 35.56 45.23                +");
   Serial.println("+          Send '0' or space to back home          +"); 
   Serial.println("+          Send 'W' to wakeup the antenna          +"); 
   Serial.println("+         Send 'S' to shutdown the antenna         +"); 
   Serial.println("+   Send 'X 52.4' will set 52.4 angle to AZ's 0    +");
   Serial.println("+   Send 'Y 30.6' will set 30.6 angle to EL's 0    +");
   Serial.println("+                 GOOD LUCK !!!                    +");
   Serial.println("+--------------------------------------------------+");
 //Serial.println("+        ____      _    ____ ___ _____ _           +");
 //Serial.println("+       |  _ \    / \  / ___|_ _| ____| |          +");
 //Serial.println("+       | |_) |  / _ \ \___ \| ||  _| | |          +");
 //Serial.println("+       |  _ <  / ___ \ ___) | || |___| |___       +");
 //Serial.println("+       |_| \_\/_/   \_\____/___|_____|_____|      +");
   Serial.println("+                  Version 2.0                     +");
   Serial.println("+--------------------------------------------------+");
   Serial.println("++++++++++++++++++++++++++++++++++++++++++++++++++++");
   
}
//////////////////////////////////////////////////////////////////////////////////////////
/*   static int antennapower(struct pt *pt) 
     {  
      PT_BEGIN(pt);  
      while(1) 
        {  
         if( power == 1 )
           {             
              digitalWrite(2,HIGH);            
           } 
         else
           {
              //delay(6000); 
              digitalWrite(2,LOW);
           }
       } 
      PT_END(pt); 
    }
*/ 
//////////////////////////////////////////////////////////////////////////////////////////

void loop()                                          //主程序循环
{

//starttime = millis();                                //计算程序循环一次耗时
//power(&pt2);
//下面为初始化各种变量，请勿移动位置
char a[10]={0};
char b[10]={0};
char c[10]={0};
char d[10]={0};
char e[10]={0};
char f[10]={0};
//char g[10]={0};
//char h[10]={0};
//char i[10]={0};
//char j[10]={0};

String tmp_a = "";
String tmp_b = "";
String tmp_c = "";
String tmp_d = "";
String tmp_e = "";
String tmp_f = "";
//String tmp_g = "";
//String tmp_h = "";
//String tmp_i = "";
//String tmp_j = "";
String usbdata = "";


//下面为接收串口数据，将串口字符串数据，传递给usbdata变量

   while (Serial.available() > 0) 
   {
    usbdata += char(Serial.read());
    delay(2);
   }                                                                    

//usb串口的数据格式例子：
//Easycomm I COM port模式的数据为：（注册后的模式）
//AZ84.46 EL19.41\#13
//UP0 DN12635455 UMFM-N DMFM-N AZ253.62 EL56.13 SNNOAA-18\#0\#13\#10
//AZ84.46 EL19.41\#13
//AZ:23.0,EL:56.0
//25.3 36.2


//自动追踪开始
 if ( usbdata.length() > 4 &&　usbdata.startsWith("A")　)
   {
      const char *usb=usbdata.c_str();                               //将string数据转成字符串数据
      //Serial.println(usb);

      //下面为追踪协议的选择，根据使用控制软件不同修改
      sscanf(usb, "AZ%s EL%s\#13",&a, &b);                           //WXtrack 的Easycomm协议
      //sscanf(usb, "AZ:%[^','],EL:%s",&a, &b);                      // split out AZ and EL yes---DDE client 
      //sscanf(usb, "UP0 DN%s UMFM-N DMFM-N AZ%s EL%s SNNOAA-%s\#0\#13\#10",&g,&a,&b,&h);    //no split  ---DDE client  
      tmp_a = a;                                                
      tmp_b = b;                                                 
                                          
      X = tmp_a.toFloat();                                            //进行浮点数转换
      Y = tmp_b.toFloat();                                      
      //Serial.println("AZ%.2f EL%.2f",X ,Y); 
      
      if( X - angle_x_tmp > 359 && angle_x == 0 )                     //检测是否是从0度过渡到359度，不至于方位角倒转一圈
         {
             angle_x = -1;
          } 
      else if( angle_x_tmp - X > 2 && angle_x == -1 )                 //检测是否是WXtrack的park antenna命令,0--> 359.x
         {
             angle_x = 0;      
          }
      else if( angle_x_tmp - X > 2 && angle_x == 1 )                  //检测是否是WXtrack的park antenna命令,359.x--> 0 
         {  
             angle_x = 0;
          }
      else if( X - angle_x_tmp < -359 && angle_x == 0 )               //检测是否是从359度过渡到0度，不至于方位角倒转一圈
         {
             angle_x = 1;           
          }

      angle_x_tmp = X;                                                //重新赋值新的比较角度值
      gotoangle_x = X;
      gotoangle_y = Y;
      _angle_1 = true;   
      _angle_2 = true;   
        
   }
//程序自动追踪代码结束
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  else if (usbdata.startsWith("W") && power == 0 && power_tmp == 1 )                                 //检测是否为天线唤醒命令“W”
  {     
      Serial.println("Antenna's power is ON, Now please wait 2 senconds !!");
      Serial.println("----------------------------------------------------"); 
      digitalWrite(power1Pin, HIGH);
      digitalWrite(power2Pin, HIGH);
      power = 1; 
      power_tmp = 0;
      delay(2000);
      Serial.println("AZ= 0.00");
      Serial.println("EL= 0.00");
      Serial.println("-----------");
      //gotoangle_x = 0;
      //gotoangle_y = 0;
      //_angle_1 = true;
      //_angle_2 = true;
      //angle_x = 0;

  }
 else if ( usbdata.startsWith("U") && power == 0 && power_tmp == -1 )                                 //检测是否为天线解除锁定命令“W”
  {     
      Serial.println("Antenna is UNLOACKED, Now please wait 2 senconds !!");
      Serial.println("---------------------------------------------------"); 
      digitalWrite(power1Pin, HIGH);
      //digitalWrite(power2Pin, HIGH);
      power = 1; 
      power_tmp = 0;
      delay(2000);
      Serial.println("OK !");
      Serial.println("-----------");

  }
 else if ( usbdata.startsWith("S") )                                                                   //检测是否为天线停止命令“S”   
  {
      power_tmp = 1; 
      gotoangle_x = 0;
      gotoangle_y = 0;
      _angle_1 = true;
      _angle_2 = true;
      angle_x = 0;  
  } 
 else if ( usbdata.startsWith("L") )                                                                   //检测是否为天线锁定命令“L”
  {
      power_tmp = -1;         
  } 
 else if ( usbdata.length() > 0 && power_tmp == 1 )
  {
         Serial.println("#############################################################"); 
         Serial.println("The antenna's power is OFF , Please send 'W' to wake him up !");
         Serial.println("#############################################################"); 
  }  
 else if ( usbdata.length() > 0 && power_tmp == -1 )   
  {
         Serial.println("######################################################"); 
         Serial.println("The antenna is LOCKED , Please send 'U' to UNLOCK him !");
         Serial.println("######################################################"); 
  }
 
 else if ( usbdata.length() > 0 && power == 1 && power_tmp == 0 )
  {                                                              
    if ( usbdata.startsWith("X") )                     //检测是否为校准方位角命令
      { 
      const char *usb=usbdata.c_str(); 
      sscanf(usb, "X%s",&e);
      tmp_e = e;
      N = tmp_e.toFloat()*AZFACTOR;
      stepper1.setCurrentPosition(N);
      _angle_1 = true;
      //_angle_2 = true;  
      }
    else if ( usbdata.startsWith("Y") )                //检测是否为校准仰角命令
      { 
      const char *usb=usbdata.c_str();
      sscanf(usb, "Y%s",&f);
      tmp_f = f;
      M = tmp_f.toFloat()*ELFACTOR;
      stepper2.setCurrentPosition(-M);
      //_angle_1 = true;
      _angle_2 = true;
      }
    //否则进行手动控制
    else                             
      {
      const char *usb=usbdata.c_str();                           //将string数据提取仰角数据
      sscanf(usb, "%s %s",&c, &d);
      tmp_c = c;
      tmp_d = d;
      X = tmp_c.toFloat();
      Y = tmp_d.toFloat();
      if (X < 400 && X > -450)                                   //手动控制时，进行角度限制。
       {        
          gotoangle_x = X;
          Serial.print("AZ= %.2f", X);
       }
      else
       {
          Serial.println("Wrong angle of AZ !");
       }
      if (Y > -20 && Y < 100)       
       {
          gotoangle_y = Y;
          Serial.print("EL= %.2f", Y);
       }
      else
       {
          Serial.println("Wrong angle of EL !");
       }
      Serial.println("-----------");
      angle_x = 0;      
      _angle_1 = true;
      _angle_2 = true;
     }
  }
                                        
//方位角旋转角度转换

   if (_angle_1)     
     {
       if ( angle_x == -1 )                              //判断是否是从0度过渡到359.x度
        {       
        stepper1.moveTo((360+(-gotoangle_x)) * AZFACTOR);   
        }
       else if ( angle_x == 1 )                          //判断方位角从359.x过渡到0度   
        {    
         stepper1.moveTo(((-gotoangle_x)-360) * AZFACTOR);          
        }
       else if ( angle_x == 0 )                      
        {       
        stepper1.moveTo((-gotoangle_x) * AZFACTOR);   
        }                                                 
     }
//仰角旋转角度转换
   
   if (_angle_2)   
       stepper2.moveTo(gotoangle_y * ELFACTOR);   
       stepper1.run();
       stepper2.run();

       

//stoptime = millis();                                              //读取程序循环结束时间参数   

if ( power_tmp == 1 )                                             //执行天线断电程序
  {
    //looptime = stoptime - starttime;
    accumulate += looptime;
    //Serial.println(accumulate);                                 //打印程序运行一次需要的时间，实测50多微秒                     
    if( accumulate > 60000 )
        {
          //Serial.println(accumulate);
          digitalWrite(power1Pin, LOW);
          digitalWrite(power2Pin, LOW);
          Serial.println("Antenna's power is OFF!");
          Serial.println("#######################");
          accumulate = 0;                                        //重置时间计数器
          power = 0; 
        }
  }

if ( power_tmp == -1 && power == 1 )
  {     
      digitalWrite(power1Pin, LOW);
      Serial.println("Antenna is LOCKED!");
      Serial.println("##################");
      power = 0;
  }

if ( power_tmp == 0 && power == 0)
   {
      digitalWrite(power1Pin, HIGH);
      digitalWrite(power2Pin, HIGH);
      power = 1;  
      delay (2000);      
   }
}
//loop循环函数结束
/*
if ( power == 1 && power_tmp == 0)
   {
      digitalWrite(power1Pin, LOW);
      delay (2000); 
   }
*/
//antennapower(&pt1);                                               //执行线程1
//xxxxxxxxxxxx(&pt2);                                               //执行线程1
  



/* 
static int protothread1(struct pt *pt) 
{  
  PT_BEGIN(pt);  
  while(power == 0) 
  {  
    delay(6000); 
    digitalWrite(12,LOW); 
  } 
  PT_END(pt); 
}
//外部中断程序
/*
void （power1）
  {
    delay(6000);
    digitalWrite(power1Pin, LOW);    
  }
  void (power2)

  {
    
  }
//整个程序结束
*/
