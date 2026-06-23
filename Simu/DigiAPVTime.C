// ------------------------------------------------------------------------------------------------------- //
// This helper script will generate GEM APV timing offsets for digitized GEM ADC samples. It is modified   //
// from Anu's script which converted gem pedestal and common mode DAQ files into a more usable format for  //
// the digitizer for added CM fluctations in simulation.                                                   //
//                                                                                                         //
// Some references to this original script (i.e. labeling of "cm") still exist.                            //
//                                                                                                         //
// Usage:                                                                                                  //
//   root -l                                                                                               //
//   .L DigiAPVTime.C+                                                                                     //
//   DigiAPVTime("db_local.txt", "APV_time.txt", "output.txt");                                            //
//                                                                                                         //
// Arguments:                                                                                              //
//   db_local : File containing APV mapping as exists in DB files for either FT or FPP                     //
//   APV_time : Output of CompareAPV.C, text file containing APV time offsets for either FT or FPP         //
//   output   : Output file  for storing results                                                           //
//                                                                                                         //
// ---------                                                                                               //
//  Jacob McMurtry, rby2vw@virginia.edu CREATED 6-17-2026                                                  //
// ---------                                                                                               //
// ** Do not tamper with this sticker! Log any updates to the script above.                                //
// ------------------------------------------------------------------------------------------------------- //


#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include "TString.h"

using namespace std;

// Parameters to be read-in from the local database.
int db_nmodules {-1};
std::vector<int> db_modAPVmap {-1};

//module id, axis, position on axis
struct apvInfoGEM {
	int gemid, axis, pos, invert;

	void print() const {
	std::cout 
	<< "GEMId: " << gemid 
	<< " Axis: " << axis 
	<< " Pos: " << pos 
	<< " Invert: " << invert
	<< std::endl;
	};
	bool operator<(const apvInfoGEM& other) const {
	if (gemid != other.gemid) return gemid < other.gemid;
	if (axis != other.axis) return axis < other.axis;
	return pos < other.pos;
	}
};

//vtpcrate, fiber, adc_ch
struct apvInfoDAQ{
	//fill these values they id the APV
	int vtpcrate;
	int fiber; //NOTE: == mpd_id
	int adc_ch;
	// int invert;

	bool operator<(const apvInfoDAQ& other) const {
        if (vtpcrate != other.vtpcrate) return vtpcrate < other.vtpcrate;
        if (fiber != other.fiber) return fiber < other.fiber;
        return adc_ch < other.adc_ch;
    }

    void print() const {
        std::cout << "VTPcrate: " << vtpcrate 
                  << " Fiber: " << fiber 
                  << " ADC_ch: " << adc_ch 
                  << std::endl; 
	}
};

// Map with each APV's DAQ info (vtpcrate, fiber, adc_ch) being used as keys to get their values (gemid, axis, pos, invert) at the GEM module's end.
std::map<apvInfoDAQ, apvInfoGEM> db_apvInfoMap;

struct CommonMode {
  double mean{};
  double sigma{};
};

using cmmap = std::map<int, CommonMode>; // Map CM mean and sigma values to each APV number.
//using trackercmmap = std::map<int, cmmap>; // Map each APV to GEM planes. GEM plane# = GEM_ID*2 + Axis_no


int readlocalDB ( const string& db_name );
int makeCMfile ( const string& cmfile_name, const string& cmfileout_name );


int DigiAPVTime( const string db_name, const string cmfile_name, const string digicmoutfile_name){

	if ( readlocalDB(db_name) < 0 ) return -1;
    if ( makeCMfile(cmfile_name,digicmoutfile_name) < 0 ) return -1;

    return 0;
}



int readlocalDB(const std::string& db_name) {
    std::cout << "Reading Local Database file " << db_name << std::endl;
    std::ifstream m_dbFile(db_name);
    if (!m_dbFile.is_open()) {
        std::cerr << "ERROR: Could not open the DB file " << db_name << std::endl;
        return -1;
    }

    TString currentline;
    // First pass: find nmodules
    while (currentline.ReadLine(m_dbFile)) {
        if (!currentline.BeginsWith("#")) {
            TObjArray* tokens = currentline.Tokenize(" ");
            if (tokens->GetEntries() > 1) {
                TString skey = ((TObjString*)(*tokens)[0])->GetString();
                if (skey == "nmodules") {
                    db_nmodules = ((TObjString*)(*tokens)[1])->GetString().Atoi();
                }
            }
        }
    }

    if (db_nmodules <= 0) {
        std::cerr << "The number of modules not defined or set to 0 in the DB " << db_name << "!!!" << std::endl;
        return -1;
    }

    db_modAPVmap.clear();
    db_modAPVmap.resize(db_nmodules);

    m_dbFile.clear();
    m_dbFile.seekg(0); // Reset file pointer for second pass

    while (currentline.ReadLine(m_dbFile)) {
        if (!currentline.BeginsWith("#")) {
            TObjArray* tokens = currentline.Tokenize(" ");
            if (tokens->GetEntries() >= 1) {
                TString skey = ((TObjString*)(*tokens)[0])->GetString();
                for (int i = 0; i < db_nmodules; i++) {
                    if (skey == Form("m%i.chanmap", i)) {
                        while (currentline.ReadLine(m_dbFile)) {
                            if (currentline.BeginsWith("#")) break;
                            TObjArray* tokens_chanmap = currentline.Tokenize(" ");
                            if (tokens_chanmap->GetEntries() == 9 &&
                                ((TObjString*)(*tokens_chanmap)[0])->GetString().IsDigit()) {
                                int vtpcrate = ((TObjString*)(*tokens_chanmap)[0])->GetString().Atoi();
                                int fiber    = ((TObjString*)(*tokens_chanmap)[2])->GetString().Atoi();
                                int adc_ch   = ((TObjString*)(*tokens_chanmap)[4])->GetString().Atoi();
                                int gemid    = i;
                                int pos      = ((TObjString*)(*tokens_chanmap)[6])->GetString().Atoi();
                                int invert   = ((TObjString*)(*tokens_chanmap)[7])->GetString().Atoi();
                                int axis     = ((TObjString*)(*tokens_chanmap)[8])->GetString().Atoi();

                                apvInfoGEM thisAPVgeminfo{gemid, axis, pos, invert};
                                apvInfoDAQ thisAPVdaqinfo{vtpcrate, fiber, adc_ch};
                                db_apvInfoMap[thisAPVdaqinfo] = thisAPVgeminfo;
                            } else break;
                            delete tokens_chanmap;
                        }
                    }                    
                }
            }
            if (tokens->GetEntries()>1){
               TString skey = ((TObjString*)(*tokens)[0])->GetString(); 
               for (int i=0; i<db_nmodules; i++){
                if(skey == Form("m%i.apvmap",i)){
                    int thismodAPVmap = ((TObjString*)(*tokens)[1])->GetString().Atoi();
                    db_modAPVmap[i] = thismodAPVmap;
                }
               }
            }
            delete tokens;
        }
    }    
    return 0;
}


int makeCMfile(const std::string& cmfile_name, const string& cmfileout_name) {
    std::cout << "Reading CM file " << cmfile_name << std::endl;
    std::ifstream m_cmfile(cmfile_name);
    if (!m_cmfile.is_open()) {
        std::cerr << "ERROR: could not open CM file " << cmfile_name << std::endl;
        return -1;
    }

    std::map<int, cmmap> trackercmmap; // digiplane -> (pos -> CM info)
    TString currentline;

    while (currentline.ReadLine(m_cmfile)) {
        if (!currentline.BeginsWith("#")) {
            TObjArray* tokens = currentline.Tokenize(" ");
            if (tokens->GetEntries() >= 5 && ((TObjString*)(*tokens)[0])->GetString().IsDigit()) {
                int vtpcrate = ((TObjString*)(*tokens)[1])->GetString().Atoi();
                int fiber    = ((TObjString*)(*tokens)[2])->GetString().Atoi();
                int adc_ch   = ((TObjString*)(*tokens)[3])->GetString().Atoi();
                double cm_mean  = ((TObjString*)(*tokens)[4])->GetString().Atof();
                double cm_sigma = ((TObjString*)(*tokens)[4])->GetString().Atof();
                
                // Lookup GEM info
                apvInfoDAQ thisAPVdaqinfo{vtpcrate, fiber, adc_ch};
                auto it = db_apvInfoMap.find(thisAPVdaqinfo);
                if (it == db_apvInfoMap.end()) {
                    std::cerr << "Warning: APV not found in DB for crate=" << vtpcrate
                              << " fiber=" << fiber << " adc_ch=" << adc_ch << std::endl;
                    delete tokens;
                    continue;
                }

                int gemid = it->second.gemid;
                int axis  = it->second.axis;
                int pos   = it->second.pos;
                
                int digiplane_num = gemid * 2 + axis;
                trackercmmap[digiplane_num][pos] = CommonMode{cm_mean, cm_sigma};
            }
            delete tokens;
        }
    }

    const double DEFAULT_CM = 0.0;
    const double DEFAULT_SIGMA = 999.0;

    for (const auto& [daqinfo, geminfo] : db_apvInfoMap) {

        int digiplane_num = geminfo.gemid * 2 + geminfo.axis;
        int pos = geminfo.pos;

        auto& apvmap = trackercmmap[digiplane_num];

        if (apvmap.find(pos) == apvmap.end()) {
            apvmap[pos] = CommonMode{DEFAULT_CM, DEFAULT_SIGMA};
        }
    }

    std::ofstream m_cmdigifile(cmfileout_name);
    m_cmdigifile << "#For each GEM_plane#:\n";
    m_cmdigifile << "# APV# time_offset\n" << '\n';

    for (const auto& [digiplane_num, apvcmmap] : trackercmmap) {
        m_cmdigifile << "GEM_plane# " << digiplane_num << "\n";
        for (const auto& [apvnum, cminfo] : apvcmmap) {
            m_cmdigifile << apvnum << " " << cminfo.mean << "\n";
        }
        m_cmdigifile << '\n';
    }

    m_cmdigifile.close();
    std::cout << "Output digi CM file " << cmfileout_name << " created.\n";
    return 0;
}