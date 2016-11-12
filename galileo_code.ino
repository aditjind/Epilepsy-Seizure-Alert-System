#include<SD.h>
#include<SPI.h>
File myFile;
int sensorPin = A0;
double seizure_data[3000]={0};
int ledPin = 13;
int sensorValue = 0;
int val;
int flag_sms=0;
// weights in linear kernel
float weights[]=
{ 20.8061,
 -221.703,
  14.4422,
  177.097,
 -406.268,
  237.769,
  -18.844,
 -566.564,
  42.9924,
 -262.076,
 -261.008,
  906.129,
  111.471,
  53.7562,
  0.343415,
  7.29209,
 -1327.42,
  1614.02,
 -486.806,
 -357.666
 };
int toggle=1,button=8,out=6,out2=7,out3=4;
void setup() 
{
    myFile = SD.open("Data_10_12_2014.txt", FILE_WRITE);
    Serial.begin(9600);
    while(!Serial)
    {
    }
//    Serial.println("Done!!!\n");
    while(!SD.exists("Data_10_12_2014.txt"))
    {
    }
    pinMode(ledPin, OUTPUT);
    pinMode(button, INPUT);
    pinMode(out, OUTPUT);
    pinMode(out2, OUTPUT);    
    digitalWrite(out,LOW);
    digitalWrite(out2,HIGH);
//    Serial.println("sdf\n");
}

void loop()
{ 
//  Serial.print(digitalRead(button));
  if(digitalRead(button)==0)
  {
//    Serial.println("enter\n");
    digitalWrite(ledPin,LOW);
    toggle=1;
    analogReadResolution(12);
    val=analogRead(sensorPin);
//    Serial.print("digitalRead(button)");
    Serial.println(val);
    myFile.println(val);
//    Serial.println("yeah");
    digitalWrite(ledPin,HIGH);
//    delay(100);
  }
  else
  {  
//      Serial.print("Yeah\n");
      if(toggle)
      {
        Serial.print("blah\n");
        myFile.close();
        myFile = SD.open("testData.txt", FILE_WRITE);
        while(!SD.exists("testData.txt"))
        {
          Serial.print("LOL");
        }
        toggle=0;
        double eeg_vect[255];
        seizure_data_plot(myFile);
        int j_val=255,index=0,i;
        for(i=0;i<=2560;++i)
        {
          eeg_vect[index]=seizure_data[i];
          if(i==j_val)
          {
             feature(eeg_vect);
             j_val+=256;
             index=-1;
             Serial.println("index is \n");
             Serial.println(-1);
          }
          ++index;
       //   Serial.println("index is \n");
       //   Serial.println(index);
        }
      }
  }
}
void seizure_data_plot(File myFile)
{
    int k=1;
    int i=0;
    double n='0';
    int c;
    int count=0;
    while (myFile.available()) 
    {
        c=myFile.read();
        while(c==' ' || c=='  ' || c=='\n')
        {
          c=myFile.read();
        }
        if(c==',')
        {
          n='0';
          k=1;
        }
        else if(c=='-')
        {
          k=-1;
        }
        else if(c=='.')
        {
          n+=(myFile.read()-48)*0.1+(myFile.read()-48)*0.01;
          n*=k;
          Serial.println(n);
          seizure_data[i++]=n;
        }
        else
        {
          if(n!='0')
          {
            count++; 
            n=n*10+(c-48);
          }
          else
          {
            n=c-48;
          }
        }
    }
    myFile.close();
}
void feature(double eeg[])
{
    //sampling rate 256 hz
    //1 sec eeg recording
    double a2[65],a3[33],a4[17],d3[33],d4[17];    //approximation and detail wavelet coefficienrs
    int j,i;
    double mean_d3=0,mean_d4=0,mean_a4=0;
    double var_d3=0,var_d4=0,var_a4=0;
    double energy_d3=0,energy_d4=0,energy_a4=0;
    double max_d3=-1000000000,max_d4=-1000000000,max_a4=-1000000000;
    double min_d3=1000000000,min_d4=1000000000,min_a4=1000000000;
    double entropy_d3=0,entropy_d4=0,entropy_a4=0;
    int len=256;    //sampling rate of ADC
    j=1;
    for(i=1;i<=len;i=i+4)
    {
         a2[j++]  =(eeg[i]+eeg[i+1]+eeg[i+2]+eeg[i+3])*0.5;  
    }
    len=64;
    double mag1,mag2;
    j=1;
    i=1;
    for(i=1;i<=len;i=i+2)
    {
         a3[j]=( a2[i] + a2[i+1] )*0.707106781;
         d3[j]=( a2[i] - a2[i+1] )*0.707106781;
         //Serial.println(a3[j]);
         if(min_d3>d3[j])
         {
             min_d3=d3[j];
         }
         if(max_d3<d3[j])                  
         {
             max_d3=d3[j];
         }
         mag1=d3[j]*d3[j];
         energy_d3+=mag1;
         entropy_d3+=log(0.0001+mag1)*mag1;
         mean_d3+=d3[j];
         j++;
    } 
    var_d3  =(energy_d3-mean_d3)/32;
    mean_d3/=32;
    len=32;
    j=1;
    for(i=1;i<=len;i=i+2)
    {
         a4[j]=(a3[i]+a3[i+1])/sqrt(2);
         d4[j]=(a3[i]-a3[i+1])/sqrt(2);
         if(min_a4>a4[j])
         {
            min_a4=a4[j];
         }
         if(max_a4<a4[j])                  
         {
             max_a4=a4[j];
         }
         if(min_d4>d4[j])
         {
           min_d4=d4[j];
         }
         if(max_d4<d4[j])
         {
             max_d4=d4[j];
         }
         mag1=d4[j]*d4[j];
         mag2=a4[j]*a4[j];         
         energy_d4+=mag1;
         energy_a4+=mag2;
         entropy_d4+=log(0.0001+mag1)*mag1;
         entropy_a4+=log(0.0001+mag2)*mag2;
         mean_d4+=d4[j];
         mean_a4+=a4[j];
         j++;
    }
    var_d4=(energy_d4-mean_d4)/16;
    var_a4=(energy_a4-mean_a4)/16;
    mean_d4/=16;    
    mean_a4/=16;
    
  
    //    calculation for mean absolute deviation
    int mad=0,iqr=0,sum=0,temp;
    for(i=1;i<=256;i++)
    {
      sum+=eeg[i];
    }
    sum/=256;
    for(i=1;i<=256;i++)
    {
      temp=eeg[i]-sum;
      if(temp<0)
      {
        temp*=-1;
      }
      mad+=temp;
    }
    //calculation for IQR(75th percentile-25th percentile)
    double var;
        //array sort
    for(i=1;i<=256;i++)
    {
        for(j=i;j<=256;j++) 
        {
          if(eeg[i]>eeg[j])
          { 
            //swap them
            var=eeg[i];
            eeg[i]=eeg[j];
            eeg[j]=var;
          }
        }
    }
    mad/=256;
    double eeg_feature[]= {mean_d3,var_d3,energy_d3,entropy_d3,max_d3,min_d3,
                           mean_a4,var_a4,energy_a4,entropy_a4,max_a4,min_a4,
                           mean_d4,var_d4,energy_d4,entropy_d4,max_d4,min_d4,
                           (eeg[192]-eeg[63]),mad
                          };
//    Serial.println("Ser vect");
//    for(i=0;i<20;i++)
//    {
//      Serial.println(eeg_feature[i]);
//    }
//    Serial.println("Ser vect end");
    double temp2=0;
    for(i=0;i<20;i++)
    {
      temp2+=weights[i]*eeg_feature[i];
    }
    if(temp2-1768921.73)
    {
      if(!flag_sms)
      {
        flag_sms=1;
        pinMode(4,OUTPUT);
        digitalWrite(4, HIGH);
        Serial.print("dsf\n");
      }
    }
}
