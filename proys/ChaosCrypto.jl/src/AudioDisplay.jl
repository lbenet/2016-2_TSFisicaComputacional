module AudioDisplay

export audioplayer

using WAV
 
function audioplayer(filepath)
  markup = """<audio controls="controls" {autoplay}>
              <source src="$filepath" />
              Your browser does not support the audio element.
              </audio>"""
 
  display(MIME("text/html") ,markup)
end
 
function audioplayer(s, fs)
    buf = IOBuffer()
    wavwrite(s, buf; Fs=fs)
    data = base64(bytestring(buf)) 
    markup = """<audio controls="controls" {autoplay}>
                <source src="data:audio/wav;base64,$data" 
                type="audio/wav" />
                Your browser does not support the audio element.
                </audio>"""
    display(MIME("text/html") ,markup)
end

end
