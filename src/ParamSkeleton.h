#include <nlohmann/json.hpp>

const nlohmann::json skeleton = R"(
{
  "mri": {
    "ute1name": "UTE",
    "ute1TR": 0.1,
    "ute2name": "UTE",
    "ute2TR": 0.1,    
    "mrac": "MRAC"
  },
  "dest": {
    "type": "VOLUME",
    "filetype": ".nii.gz",
    "path": "."
  },
  "logpath" : ".",
  "reg": {
    "name": "ANTS",
    "args": "cli args",
    "templatepath": "",
    "maskpath": ""
  }
}
)"_json;